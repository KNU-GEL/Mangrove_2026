#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import json
import requests
from requests.auth import HTTPBasicAuth
import pathlib
import time
import errno
import multiprocessing
from datetime import datetime, timedelta
import re
import random


# ============================================================
# USER CONFIGURATION
# ============================================================

# --- API ---
PLANET_API_KEY = "PLANET_API_KEY"
BASE_URL       = "https://api.planet.com/data/v1"
ORDERS_URL     = "https://api.planet.com/compute/ops/orders/v2"

# --- Paths ---
GEOMETRY_PATH  = "/home/hjlee/project/mangrove/geojson"
OUTPUT_DIR     = "/home/hjlee/project/mangrove/data"

# --- Time Range ---
MIN_YEAR       = 2018
MAX_YEAR       = 2025

# --- Parallel ---
NUM_CORES      = 10

# --- Imagery Options ---
N_BANDS        = 4          # 4 or 8
HARMONIZE      = True       # True: Sentinel-2 / False: no harmonization
CLOUD_PCT      = 0.8        # 0.0 ~ 1.0

# --- Order Settings ---
CHUNK_SIZE     = 200
ORDER_TYPE     = "partial"  # "full" or "partial"
CHECK_EXISTING_ORDERS = True

# ============================================================
# DERIVED SETTINGS (auto-configured from N_BANDS)
# ============================================================

if N_BANDS == 8:
    ASSET_FILTER    = ["ortho_analytic_8b_sr", "ortho_udm2"]
    PRODUCT_BUNDLE  = "analytic_8b_sr_udm2"
    INSTRUMENT_FILTER = ["PSB.SD"]
else:
    ASSET_FILTER    = ["ortho_analytic_4b_sr", "ortho_udm2"]
    PRODUCT_BUNDLE  = "analytic_sr_udm2"
    INSTRUMENT_FILTER = None  # All instruments

# ============================================================

# Derived URLs
STATS_URL = "{}/stats".format(BASE_URL)
QUICK_URL = "{}/quick-search".format(BASE_URL)


def setup_filter(coords, minyear, maxyear):
    
    geometry_filter = {"type": "GeometryFilter",
        "field_name": "geometry",
        "config": {
          "type": "Polygon",
          "coordinates": coords
        }
    }
    
    date_filter = {
        "type": "DateRangeFilter",
        "field_name": "acquired",
        "config": {
            "gte": "{}-01-01T00:00:00.000Z".format(minyear),
            "lt":"{}-12-31T23:59:59.000Z".format(maxyear)
            }
        }
    

    ground_control =  {
                    "type": "StringInFilter",
                    "config": ["true"],
                    "field_name": "ground_control"
                  }
    quality_category = {
                    "type": "StringInFilter",
                    "config": ["standard"],
                    "field_name": "quality_category"
                  }
    cloud_cover =  {
                     "type": "RangeFilter",
                     "field_name": "cloud_cover",
                     "config": {"gte": 0.0, "lte": CLOUD_PCT}
                  }
    asset = {
            "type": "AssetFilter",
            "config": ASSET_FILTER
        }

    permission = {
        "type":"PermissionFilter",
        "config":[
            "assets:download"
            ]
        }

    filter_config = [ 
        cloud_cover,
        quality_category,
        ground_control,
        date_filter,
        geometry_filter,
        permission,
        asset
    ]

    # 8-band: restrict to SuperDove (PSB.SD) only
    if INSTRUMENT_FILTER is not None:
        instrument = {
            "type": "StringInFilter",
            "field_name": "instrument",
            "config": INSTRUMENT_FILTER
        }
        filter_config.append(instrument)

    and_filter = {
        "type": "AndFilter",
        "config": filter_config
        }
    
    return and_filter

# [PATCH 5] Generic retry wrapper for session requests (429/5xx/ConnectionError)
def session_request_with_retry(session, method, url, max_retries=5, **kwargs):
    """Retry session request on 429, 5xx, and connection errors. Returns response or None."""
    if 'timeout' not in kwargs:
        kwargs['timeout'] = 60
    for attempt in range(1, max_retries + 1):
        try:
            if method == 'post':
                res = session.post(url, **kwargs)
            else:
                res = session.get(url, **kwargs)
        except requests.RequestException as e:
            print("Request error (attempt {}/{}): {}".format(attempt, max_retries, e))
            if attempt < max_retries:
                time.sleep(30 * attempt)
            continue
        if res.status_code == 429:
            wait = int(res.headers.get('Retry-After', 30 * attempt))
            print("429 rate limit (attempt {}/{}), waiting {}s".format(attempt, max_retries, wait))
            time.sleep(wait)
            continue
        if res.status_code >= 500:
            print("Server error {} (attempt {}), retrying...".format(res.status_code, attempt))
            if attempt < max_retries:
                time.sleep(10 * attempt)
            continue
        return res
    print("All {} retries failed for {}".format(max_retries, url))
    return None


def read_geometry(path):
    
    if(False == os.path.exists(path)):
        sys.exit("GeoJSON path doesn't exist: {}".format(path))
        
    with open(path, "r") as file1:
        geo = json.load(file1)
        
    return geo


def safe_feature_name(x):
    """Extract feature name from properties (type-agnostic)"""
    props = x.get("properties")

    if isinstance(props, list):
        val = props[1] if len(props) > 1 else "site"
    elif isinstance(props, dict):
        vals = list(props.values())
        val = vals[1] if len(vals) > 1 else "site"
    else:
        val = props if props is not None else "site"

    s = str(val).strip()
    s = re.sub(r"\s+", "_", s)                # whitespace -> _
    s = re.sub(r"[^\w.\-]+", "_", s)          # dangerous chars -> _
    return s or "site"
    

def place_order(request, auth, order_name):
    headers = {'content-type': 'application/json'}
    
    response = requests.post(ORDERS_URL, data=json.dumps(request), auth=auth, headers=headers)
    print(response)
    print(response.reason)
    
    if(response.status_code != 202):
        print("Order failed for {}".format(order_name))
        return -1, None, response.status_code
    order_id = response.json()['id']
    print(order_id)
    order_url = ORDERS_URL + '/' + order_id
    return order_url, order_id, response.status_code


def place_order_with_retry(request, auth, order_name, max_retries=10):
    """Place order with retry on 429. Returns (order_url, order_id) or (-1, None)."""
    for attempt in range(1, max_retries + 1):
        result, order_id, status_code = place_order(request, auth, order_name)
        
        if result != -1:
            return result, order_id
        
        if status_code != 429:
            print("Non-retryable error (status {}), giving up on {}".format(status_code, order_name))
            return -1, None
        
        wait_time = 30 * attempt  # 30, 60, 90, ... 300
        print("429 rate limit. Retry {}/{} after {}s for {}".format(attempt, max_retries, wait_time, order_name))
        time.sleep(wait_time)
    
    print("Max retries reached. Giving up on {}".format(order_name))
    return -1, None


def download_results(results, output_dir, order_name, overwrite=False):
    results_urls = [r['location'] for r in results]
    results_names = [r['name'] for r in results]
    print('{} items to download'.format(len(results_urls)))
    
    failed_count = 0
    skipped_count = 0
    success_count = 0
    
    filename = "failed_downloads_{}.txt".format(order_name)
    
    with open(os.path.join(output_dir, filename), "w") as file1:
    
        for url, name in zip(results_urls, results_names):
            path = pathlib.Path(os.path.join(output_dir, 'data', name))

            if overwrite or not path.exists():
                # [PATCH 3] timeout + 429/5xx retry
                max_dl_retries = 3
                r = None
                for dl_attempt in range(1, max_dl_retries + 1):
                    try:
                        r = requests.get(url, allow_redirects=True, timeout=(30, 300), stream=True)
                    except requests.RequestException as e:
                        print('Download error for {} (attempt {}): {}'.format(name, dl_attempt, e))
                        if dl_attempt < max_dl_retries:
                            time.sleep(2 ** dl_attempt)
                        continue
                    if r.status_code == 429:
                        wait = int(r.headers.get('Retry-After', 2 ** dl_attempt * 10))
                        print('429 rate limit for {}, waiting {}s (attempt {})'.format(name, wait, dl_attempt))
                        time.sleep(wait)
                        r.close()
                        r = None
                        continue
                    if r.status_code >= 500:
                        print('Server error {} for {} (attempt {})'.format(r.status_code, name, dl_attempt))
                        r.close()
                        if dl_attempt < max_dl_retries:
                            time.sleep(2 ** dl_attempt * 5)
                        r = None
                        continue
                    break  # success or 4xx (except 429)

                if r is None or r.status_code != 200:
                    status = r.status_code if r else "NO_RESPONSE"
                    print('Status code {}, {} not downloaded.'.format(status, name))
                    if r is not None:
                        r.close()
                    failed_count += 1
                    file1.write("{}, {}, {}, {} \n".format(time.strftime("%Y%m%d-%H%M%S"), status, name, url))
                else:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    part_path = pathlib.Path(str(path) + '.part')
                    try:
                        with open(part_path, 'wb') as out:
                            for chunk in r.iter_content(chunk_size=1024*1024):
                                if chunk:
                                    out.write(chunk)
                        os.replace(part_path, path)
                        success_count += 1
                    except Exception as e:
                        print('Write error for {}: {}'.format(name, e))
                        if part_path.exists():
                            part_path.unlink()
                        failed_count += 1
                        file1.write("{}, {}, {}, {} \n".format(time.strftime("%Y%m%d-%H%M%S"), "WRITE_ERR", name, url))
                    finally:
                        r.close()
            else:
                skipped_count += 1
                
    print("\n {} Success: {} Skipped: {} Failed: {}".format(order_name, success_count, skipped_count, failed_count))
    
    return failed_count, skipped_count, success_count

def p(data):
    print(json.dumps(data, indent=2))


def rewrite_name_id_file(path, mapping):
    """Atomically rewrite a name,id file from dict. Uses tmp+replace for crash safety."""
    tmp = path + ".tmp"
    try:
        with open(tmp, "w") as f:
            for n, oid in mapping.items():
                if oid:
                    f.write("{},{}\n".format(n, oid))
        os.replace(tmp, path)
    except Exception as e:
        print("Warning: failed to rewrite {}: {}".format(path, e))
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass

# [PATCH 1] timeout - prevent infinite polling (default 2 hours)
def poll_for_success(order_url, auth, max_wait=7200):
    count = 0
    end_states = ['success', 'failed', 'partial']
    state = ''
    start = time.time()
    
    while(state not in end_states):
        # timeout check
        elapsed = time.time() - start
        if elapsed > max_wait:
            print("\n poll_for_success TIMEOUT after {:.0f}s for {}".format(elapsed, order_url))
            return None
        
        count += 1
        time.sleep(20 + random.randint(0, 5) + random.random())
        try:
            r = requests.get(order_url, auth=auth, timeout=60)
        except requests.RequestException as e:
            print("\n Poll request error: {}".format(e))
            continue
        if(r.status_code != 200):
            print("\n Status code: {}".format(r.status_code))
            continue
        
        response = r.json()
        state = response['state']
        print(state)
        if state in end_states:
            return r
            break


# Worker function for parallel processing
def process_geojson_worker(args):
    """Process assigned geojson files for each worker"""
    worker_files, worker_id, total_workers = args
    
    # Add delay for each worker to avoid API rate limit
    start_delay = (worker_id - 1) * 15  # 15 seconds between each worker
    print(f"\n[Worker {worker_id}] Starting in {start_delay} seconds...")
    time.sleep(start_delay)
    
    # Re-seed random for each worker (prevent fork seed duplication)
    random.seed()
    
    # Set environment variables for worker
    os.environ['WORKER_ID'] = str(worker_id)
    os.environ['WORKER_FILES'] = ','.join(worker_files)
    
    print(f"\n[Worker {worker_id}/{total_workers}] Processing {len(worker_files)} files")
    for f in worker_files:
        print(f"  - {f}")
    
    # Execute main function in worker mode (pass argv directly, no sys.argv mutation)
    main(['--worker'])


# Run parallel processing
# [PATCH 2] maxtasksperchild + map_async timeout - prevent Pool deadlock
def run_parallel():
    """Manage parallel processing from main process"""
    print(f"\n===== Parallel Processing Mode ({NUM_CORES} cores) =====")
    
    # Clean up previous worker files
    for f in os.listdir(OUTPUT_DIR):
        if f.startswith("completed_orders_worker") or f.startswith("placed_orders_worker") or f.startswith("partial_orders_worker") or f.startswith("permanent_partial_orders_worker"):
            os.remove(os.path.join(OUTPUT_DIR, f))
            print(f"Removed old worker file: {f}")
    
    # Get geojson file list
    paths_list = os.listdir(GEOMETRY_PATH)
    paths_list.sort()
    geojson_files = [f for f in paths_list if f.endswith(".geojson")]
    
    print(f"Total {len(geojson_files)} files to {NUM_CORES} cores")
    
    # Assign files to cores
    batches = [[] for _ in range(NUM_CORES)]
    for i, f in enumerate(geojson_files):
        batches[i % NUM_CORES].append(f)
    
    # Print assignment info
    for i, batch in enumerate(batches):
        print(f"\nWorker {i+1}: {len(batch)} files")
    
    # Execute worker processes
    # maxtasksperchild=1: worker replaced after 1 task - prevents futex deadlock
    with multiprocessing.Pool(processes=NUM_CORES, maxtasksperchild=1) as pool:
        args_list = [(batch, i+1, NUM_CORES) for i, batch in enumerate(batches) if batch]
        # map_async + timeout: 120 hour timeout : BU SCC limit
        try:
            result = pool.map_async(process_geojson_worker, args_list)
            result.get(timeout=432000)
        except multiprocessing.TimeoutError:
            print("\n===== TIMEOUT: Pool did not finish within 120 hours, terminating =====")
            pool.terminate()
            pool.join()
        except KeyboardInterrupt:
            print("\n===== KeyboardInterrupt: terminating pool =====")
            pool.terminate()
            pool.join()
            raise
    
    print("\n===== Parallel Processing Complete =====")
    
    # Merge worker files into main files (completed + placed + partial + permanent_partial)
    for file_prefix in ["completed_orders", "placed_orders", "partial_orders", "permanent_partial_orders"]:
        main_file = os.path.join(OUTPUT_DIR, "{}.txt".format(file_prefix))
        wk_files = [os.path.join(OUTPUT_DIR, f) for f in os.listdir(OUTPUT_DIR) if f.startswith("{}_worker".format(file_prefix))]
        
        if not wk_files:
            continue
        
        all_orders = {}
        
        # Read existing main file
        if os.path.exists(main_file):
            with open(main_file, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if "," in line:
                        parts = line.split(",", 1)
                        all_orders[parts[0].strip()] = parts[1].strip()
                    else:
                        all_orders[line] = ""
        
        # Read and merge worker files
        for wf in wk_files:
            with open(wf, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if "," in line:
                        parts = line.split(",", 1)
                        all_orders[parts[0].strip()] = parts[1].strip()
                    else:
                        all_orders[line] = ""
        
        # Write merged result
        with open(main_file, "w") as f:
            for name in sorted(all_orders.keys()):
                f.write("{},{}\n".format(name, all_orders[name]))
        
        # Clean up worker files
        for wf in wk_files:
            os.remove(wf)
        
        print("Merged {} {} from {} worker files".format(len(all_orders), file_prefix, len(wk_files)))


def main(argv):
    start_time = time.time()
    
    # Check parallel processing mode
    if "--parallel" in argv and "--worker" not in argv:
        # Master process: start parallel processing
        print("Running in parallel mode")
        run_parallel()
        print("--- %s seconds ---" % (time.time() - start_time))
        return
    
    # Worker process or normal execution
    print("\n" + "=" * 50)
    print("Configuration:")
    print("=" * 50)
    print("  GEOMETRY_PATH = {}".format(GEOMETRY_PATH))
    print("  OUTPUT_DIR    = {}".format(OUTPUT_DIR))
    print("  MIN_YEAR      = {}".format(MIN_YEAR))
    print("  MAX_YEAR      = {}".format(MAX_YEAR))
    print("  N_BANDS       = {}".format(N_BANDS))
    print("  HARMONIZE     = {}".format(HARMONIZE))
    print("  CLOUD_PCT     = {}".format(CLOUD_PCT))
    print("  PRODUCT_BUNDLE= {}".format(PRODUCT_BUNDLE))
    print("  ASSET_FILTER  = {}".format(ASSET_FILTER))
    print("  INSTRUMENT    = {}".format(INSTRUMENT_FILTER))
    print("  CHUNK_SIZE    = {}".format(CHUNK_SIZE))
    print("  ORDER_TYPE    = {}".format(ORDER_TYPE))
    print("=" * 50)
    
    if(os.path.isdir(OUTPUT_DIR) == False):
        sys.exit("{} does not exist.".format(OUTPUT_DIR))
    
    session = requests.Session()
    session.auth = (PLANET_API_KEY, "")
    
    # [PATCH 4] Retry initial connection (prevents worker crash on transient DNS failure)
    for _conn_attempt in range(1, 6):
        try:
            res = session.get(BASE_URL, timeout=60)
            if res.status_code == 200:
                print("Base server is alive.")
                break
        except requests.RequestException as e:
            print("Connection attempt {}/5 failed: {}".format(_conn_attempt, e))
        if _conn_attempt < 5:
            time.sleep(30 * _conn_attempt)
    else:
        print("All connection attempts failed, skipping this worker.")
        return
        
    paths_list = os.listdir(GEOMETRY_PATH)
    paths_list.sort()
    print("\n Number of sites: {}".format(len(paths_list)))
    
    geojson_files = [f for f in paths_list if f.endswith(".geojson")]
    print("Number of geojson files:", len(geojson_files))
    
    # For worker mode, process only assigned files
    if "--worker" in argv:
        worker_files_env = os.environ.get('WORKER_FILES', '').split(',')
        worker_id = os.environ.get('WORKER_ID', '?')
        geojson_files = [f for f in geojson_files if f in worker_files_env]
        print(f"\n[Worker {worker_id}] Processing {len(geojson_files)} assigned files only")
    else:
        pass
    
    # Load locally completed orders: {order_name: order_id}
    completed_orders_file = os.path.join(OUTPUT_DIR, "completed_orders.txt")
    completed_orders = {}
    
    for f_path in [completed_orders_file] + \
        [os.path.join(OUTPUT_DIR, f) for f in os.listdir(OUTPUT_DIR) if f.startswith("completed_orders_worker")]:
        if os.path.exists(f_path):
            with open(f_path, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if "," in line:
                        parts = line.split(",", 1)
                        completed_orders[parts[0].strip()] = parts[1].strip()
                    else:
                        completed_orders[line] = None
    print("Loaded {} completed orders from local records".format(len(completed_orders)))
    
    # Load locally placed orders (ordered but not yet downloaded): {order_name: order_id}
    placed_orders_file = os.path.join(OUTPUT_DIR, "placed_orders.txt")
    placed_orders = {}
    
    for f_path in [placed_orders_file] + \
        [os.path.join(OUTPUT_DIR, f) for f in os.listdir(OUTPUT_DIR) if f.startswith("placed_orders_worker")]:
        if os.path.exists(f_path):
            with open(f_path, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if "," in line:
                        parts = line.split(",", 1)
                        placed_orders[parts[0].strip()] = parts[1].strip()
    print("Loaded {} placed orders from local records".format(len(placed_orders)))
    
    # Load locally partial orders (API returned partial, not all items delivered): {order_name: order_id}
    partial_orders_file = os.path.join(OUTPUT_DIR, "partial_orders.txt")
    partial_orders = {}
    
    for f_path in [partial_orders_file] + \
        [os.path.join(OUTPUT_DIR, f) for f in os.listdir(OUTPUT_DIR) if f.startswith("partial_orders_worker")]:
        if os.path.exists(f_path):
            with open(f_path, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if "," in line:
                        parts = line.split(",", 1)
                        partial_orders[parts[0].strip()] = parts[1].strip()
    print("Loaded {} partial orders from local records".format(len(partial_orders)))
    
    # Load permanent partial orders (retry also failed, skip in future runs): {order_name: order_id}
    permanent_partial_orders_file = os.path.join(OUTPUT_DIR, "permanent_partial_orders.txt")
    permanent_partial_orders = {}
    
    for f_path in [permanent_partial_orders_file] + \
        [os.path.join(OUTPUT_DIR, f) for f in os.listdir(OUTPUT_DIR) if f.startswith("permanent_partial_orders_worker")]:
        if os.path.exists(f_path):
            with open(f_path, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if "," in line:
                        parts = line.split(",", 1)
                        permanent_partial_orders[parts[0].strip()] = parts[1].strip()
    print("Loaded {} permanent partial orders from local records".format(len(permanent_partial_orders)))
    
    # Worker mode: write to separate file
    if "--worker" in argv:
        wid = os.environ.get("WORKER_ID", "0")
        completed_orders_file = os.path.join(OUTPUT_DIR, "completed_orders_worker{}.txt".format(wid))
        placed_orders_file = os.path.join(OUTPUT_DIR, "placed_orders_worker{}.txt".format(wid))
        partial_orders_file = os.path.join(OUTPUT_DIR, "partial_orders_worker{}.txt".format(wid))
        permanent_partial_orders_file = os.path.join(OUTPUT_DIR, "permanent_partial_orders_worker{}.txt".format(wid))
    
    
    for site_num, file in enumerate(geojson_files, start=1):
        print("\n\n===== Processing file {}/{}: {} =====".format(site_num, len(geojson_files), file))
        
        file_path = os.path.join(GEOMETRY_PATH, file)
        site_name = os.path.splitext(file)[0]
        output_site_dir = os.path.join(OUTPUT_DIR, site_name) 
        
        if not os.path.exists(output_site_dir):
            try:
                print("\n Make directory")
                os.makedirs(output_site_dir)
            except OSError as exc:
                if exc.errno != errno.EEXIST:
                    raise
        
        print("Output dir updated to: {}".format(output_site_dir))
        
        geo = read_geometry(file_path)
        
        for feat_idx, x in enumerate(geo['features']):
            feature_name = safe_feature_name(x)
            
            feature_coords = x['geometry']['coordinates']
            
            filter = setup_filter(feature_coords, MIN_YEAR, MAX_YEAR)
            
            print("Filter config for search:")
            p(filter)
            
            print("\n\nRequesting stat count.")
            request = {
                "interval" : "year",
                "filter" : filter,
                "item_types" : ["PSScene"]
                }

            res = session_request_with_retry(session, 'post', STATS_URL, json=request)
            
            if(res is None or res.status_code != 200):
                print("Stats search failed, skipping this feature.")
                continue

            for bucket in res.json()['buckets']:
                print("start_time: {} count: {}".format(bucket["start_time"], bucket["count"]))
                
            print("do a real search.")
            request = {
                "filter" : filter,
                "item_types" : ["PSScene"]
                }

            res = session_request_with_retry(session, 'post', QUICK_URL, json=request)
            
            if(res is None or res.status_code != 200):
                print("Quick search failed, skipping this feature.")
                continue
                
            geojson = res.json()
            
            print("\n Creating file to save quick search results.")
            filename = "{}_quick_search_result_{}_{}.json".format(feature_name, MIN_YEAR, MAX_YEAR)
            with open(os.path.join(output_site_dir, filename), 'w') as outfile:
                outfile.write(json.dumps(geojson) + '\n')
                
            print("\n Assemble ID list for download")
            features = geojson["features"]
            
            if(len(features) == 0):
                print("0 IDs returned in quick search, skipping this AOI.")
                continue
            id_list = []
            while(len(geojson["features"]) > 0):
                
                for x in geojson["features"]:
                    id_list.append(x["id"])
                
                next_url = geojson["_links"].get("_next")
                if not next_url:
                    break 
                
                time.sleep(9 + random.randint(0, 6))
                res = session_request_with_retry(session, 'get', next_url)
                
                if(res is None or res.status_code != 200):
                    print("Next page retrieval failed, stopping pagination.")
                    break
                    
                geojson = res.json()
                with open(os.path.join(output_site_dir, filename), 'a') as outfile:
                    outfile.write(json.dumps(geojson) + '\n')
            
            print("\n Total IDs: {}".format(len(id_list)))
            
            print("\n Chunking ID list.")
            chunks = [id_list[x:x+CHUNK_SIZE] for x in range(0, len(id_list), CHUNK_SIZE)]
            
            print("\n Number of chunks: {}".format(len(chunks)))
            
            if(len(chunks) >= 80):
                sys.exit("{} Chunks which is greater than 80.  This will exceed order capacity".format(len(chunks)))
        
            print("\n Checking connection with order server...")
            auth = HTTPBasicAuth(PLANET_API_KEY, '')
            
            thirty_days_ago = (datetime.now() - timedelta(days=30)).isoformat() + "Z"
            orders_url_filtered = ORDERS_URL + "?created__gte={}".format(thirty_days_ago)
            
            response = requests.get(orders_url_filtered, auth=auth, timeout=60)
            
            if(response.status_code != 200):
                sys.exit("Failed to connect to order server with code {}".format(response.status_code))
            else:
                print("connected!")
                
            orders_list = response.json()["orders"]
            print("Found {} orders from last 30 days".format(len(orders_list)))
                
            if CHECK_EXISTING_ORDERS:
                page_count = 0
                max_pages = 10
                while("next" in response.json()["_links"] and page_count < max_pages):
                    time.sleep(9 + random.randint(0, 6))
                    next_url = response.json()["_links"]["next"]
                    response = session.get(next_url, timeout=60)
                    print("Checking page {} of existing orders...".format(page_count + 2))
                    if(response.status_code == 200):
                        if("orders" in response.json()):
                            orders_list.extend(response.json()["orders"])
                    page_count += 1
                    
                print("Total {} orders found".format(len(orders_list)))
            
            print("\n Time to place order...")
            
            orders_to_monitor = []
            bad_order_count = 0
            
            for count, chunk in enumerate(chunks):
                print("\n Processing chunk: {}".format(count))
                
                order_name = "{}_chunk_{}_{}_{}".format(site_name, count, MIN_YEAR, MAX_YEAR)
                
                # 1) completed ?? skip
                if order_name in completed_orders:
                    print("\n Already completed (downloaded), skipping: {}".format(order_name))
                    continue
                
                # 1b) permanent_partial ?? skip (retry also failed, manual review needed)
                if order_name in permanent_partial_orders:
                    print("\n Permanent partial (retry also failed), skipping: {}".format(order_name))
                    continue
                
                # 2) partial ?? missing items retry
                if order_name in partial_orders:
                    missing_file = os.path.join(output_site_dir, "missing_items_{}.json".format(order_name))
                    missing_ids = []
                    
                    if os.path.exists(missing_file):
                        with open(missing_file, "r") as f:
                            missing_ids = json.load(f)
                    else:
                        # Fallback: reconstruct missing from order.json + order_result.json
                        order_json_path = os.path.join(output_site_dir, "order_{}.json".format(order_name))
                        result_json_path = os.path.join(output_site_dir, "order_result_{}.json".format(order_name))
                        if os.path.exists(order_json_path) and os.path.exists(result_json_path):
                            print("\n Reconstructing missing items from order/result files")
                            with open(order_json_path, "r") as f:
                                order_data = json.load(f)
                            with open(result_json_path, "r") as f:
                                result_data = json.load(f)
                            requested_ids = set(order_data["products"][0]["item_ids"])
                            # FIX: substring match instead of split('/') exact match
                            delivered_ids = set()
                            for req_id in requested_ids:
                                for res_item in result_data.get('_links', {}).get('results', []):
                                    if req_id in res_item.get('name', ''):
                                        delivered_ids.add(req_id)
                                        break
                            missing_ids = list(requested_ids - delivered_ids)
                            if missing_ids:
                                with open(missing_file, "w") as f:
                                    json.dump(missing_ids, f)
                                print(" Reconstructed {} missing items".format(len(missing_ids)))
                        else:
                            print("\n Partial but cannot reconstruct missing for {}, skipping".format(order_name))
                            review_log = os.path.join(output_site_dir, "partial_needs_review.txt")
                            with open(review_log, "a") as f:
                                f.write("{},{}\n".format(time.strftime("%Y%m%d-%H%M%S"), order_name))
                            continue
                    
                    if len(missing_ids) > 0:
                        retry_order_name = "{}_retry".format(order_name)
                        if retry_order_name in completed_orders or retry_order_name in placed_orders or retry_order_name in partial_orders:
                            print("\n Retry already exists for {}, skipping".format(order_name))
                            continue
                        print("\n Partial order has {} missing items, creating retry order: {}".format(
                            len(missing_ids), retry_order_name))
                        tools = []
                        if HARMONIZE:
                            tools.append({"harmonize": {"target_sensor": "Sentinel-2"}})
                        tools.append({"clip": {"aoi": {"type": "Polygon", "coordinates": feature_coords}}})
                        
                        retry_request = {
                            "name": retry_order_name,
                            "order_type": "full",
                            "products": [{
                                "item_ids": missing_ids,
                                "item_type": "PSScene",
                                "product_bundle": PRODUCT_BUNDLE
                            }],
                            "tools": tools
                        }
                        
                        retry_filename = "order_{}.json".format(retry_order_name)
                        with open(os.path.join(output_site_dir, retry_filename), 'w') as outfile:
                            json.dump(retry_request, outfile)
                        
                        time.sleep(9 + random.randint(0, 6))
                        order_url, order_id = place_order_with_retry(retry_request, auth, retry_order_name)
                        
                        if order_url == -1:
                            bad_order_count += 1
                            print(" Retry order failed for {}".format(retry_order_name))
                            failed_log = os.path.join(output_site_dir, "failed_orders.txt")
                            with open(failed_log, "a") as f:
                                f.write("{}, RETRY_FAIL, {}, from {}\n".format(
                                    time.strftime("%Y%m%d-%H%M%S"), retry_order_name, order_name))
                        else:
                            print(" Retry order placed: {} (id: {})".format(order_url, order_id))
                            with open(placed_orders_file, "a") as f:
                                f.write("{},{}\n".format(retry_order_name, order_id))
                            placed_orders[retry_order_name] = order_id
                            orders_to_monitor.append((order_url, retry_order_name))
                            # Keep original in partial_orders until retry actually succeeds.
                            # On retry SUCCESS, original will be promoted to completed.
                    else:
                        print("\n Partial but no missing items for {}, skipping".format(order_name))
                    continue
                
                # 3) placed ?? monitor/download
                if order_name in placed_orders:
                    local_order_id = placed_orders[order_name]
                    order_url = ORDERS_URL + '/' + local_order_id
                    print("\n Found in placed record (id={}), adding to monitor".format(local_order_id))
                    orders_to_monitor.append((order_url, order_name))
                    continue
                
                # 4) Fallback: check server orders by name
                if CHECK_EXISTING_ORDERS:
                    print("\n Checking existing orders on server")
                    
                    found_existing = False
                    for order in orders_list:
                        if(order["name"] == order_name):
                            print("\n Found existing Order on server by name.")
                            order_url = order["_links"]["_self"]
                            server_order_id = order_url.rsplit('/', 1)[-1]
                            with open(placed_orders_file, "a") as f:
                                f.write("{},{}\n".format(order_name, server_order_id))
                            placed_orders[order_name] = server_order_id
                            orders_to_monitor.append((order_url, order_name))
                            found_existing = True
                            break
                    
                    if found_existing:
                        continue
                
                # 5) New order
                feature_coords_buffer = feature_coords
                
                tools = []
                if HARMONIZE:
                    tools.append({
                        "harmonize": {
                            "target_sensor": "Sentinel-2"
                        }
                    })
                tools.append({
                    "clip": {
                        "aoi": {
                            "type": "Polygon",
                            "coordinates": feature_coords_buffer
                        }
                    }
                })
                
                request = {  
                   "name": order_name,
                   "order_type": ORDER_TYPE,
                   "products":[
                      {  
                         "item_ids": chunk,
                         "item_type": "PSScene",
                         "product_bundle": PRODUCT_BUNDLE
                      }
                   ],
                   "tools": tools
                }
                
                filename = "order_{}.json".format(order_name)
                with open(os.path.join(output_site_dir, filename), 'w') as outfile:
                    json.dump(request, outfile)
                
                time.sleep(9 + random.randint(0, 6))
                
                order_url, order_id = place_order_with_retry(request, auth, order_name)
                
                if(order_url == -1):
                    bad_order_count = bad_order_count +1
                    failed_log = os.path.join(output_site_dir, "failed_orders.txt")
                    with open(failed_log, "a") as f:
                        f.write("{}, {}\n".format(time.strftime("%Y%m%d-%H%M%S"), order_name))
                    continue
                else:
                    print("\n Placed order. Order url: {} (id: {})".format(order_url, order_id))
                    with open(placed_orders_file, "a") as f:
                        f.write("{},{}\n".format(order_name, order_id))
                    placed_orders[order_name] = order_id
                    orders_to_monitor.append((order_url, order_name))
            
            print("\n {} out of {} orders not placed successfully.".format(bad_order_count, len(chunks)))      
            
            if(len(orders_to_monitor) == 0):
                print("No orders to monitor for {}, skipping.".format(site_name))
                continue
            
            total_files_to_download = 0
            total_failed = 0
            total_skipped = 0
            total_success = 0
            order_list_failed = []
            order_list_failed_download = []
            
            print("\n Monitor orders...")
            for num, (order_url, order_name) in enumerate(orders_to_monitor):
                try:
                    print("Monioring order num {}, name: {}, url: {}".format(num, order_name, order_url))
                 
                    r = poll_for_success(order_url, auth)
                    
                    if r is None:
                        print("\n Order poll timed out for {}, skipping".format(order_name))
                        order_list_failed.append(num)
                        continue
                    
                    response = r.json()
                    state = response['state']
                    
                    # FAILED
                    if state == "failed":
                        print("\n Order FAILED, removing from placed to allow re-order: {}".format(order_name))
                        placed_orders.pop(order_name, None)
                        rewrite_name_id_file(placed_orders_file, placed_orders)
                        order_list_failed.append(num)
                        continue
                    
                    # PARTIAL
                    if state == "partial":
                        print("\n Order PARTIAL: not all items delivered by Planet.")
                        results = response['_links']['results']
                        print(" Delivered files: {}".format(len(results)))
                        total_files_to_download += len(results)
                        
                        filename = "order_result_{}.json".format(order_name)
                        with open(os.path.join(output_site_dir, filename), 'w') as outfile:
                            json.dump(response, outfile)
                        
                        failed_count, skipped_count, success_count = download_results(results, output_site_dir, order_name)
                        
                        # Identify missing items: requested vs delivered (substring match)
                        order_json_path = os.path.join(output_site_dir, "order_{}.json".format(order_name))
                        missing_ids = []
                        if os.path.exists(order_json_path):
                            with open(order_json_path, "r") as f:
                                order_data = json.load(f)
                            requested_ids = set(order_data["products"][0]["item_ids"])
                            # FIX: substring match instead of split('/') exact match
                            delivered_ids = set()
                            for req_id in requested_ids:
                                for res_item in results:
                                    if req_id in res_item.get('name', ''):
                                        delivered_ids.add(req_id)
                                        break
                            missing_ids = list(requested_ids - delivered_ids)
                            print(" Requested: {}, Delivered: {}, Missing: {}".format(
                                len(requested_ids), len(delivered_ids), len(missing_ids)))
                        
                        if missing_ids:
                            missing_file = os.path.join(output_site_dir, "missing_items_{}.json".format(order_name))
                            with open(missing_file, "w") as f:
                                json.dump(missing_ids, f)
                            print(" Missing items saved to: {}".format(missing_file))
                        
                        order_id = placed_orders.get(order_name) or order_url.rsplit("/", 1)[-1]
                        
                        # If this is a retry order that is also partial, move to permanent_partial
                        if order_name.endswith("_retry"):
                            original_name = order_name[:-len("_retry")]
                            print("\n Retry also PARTIAL. Giving up on {} and original {}".format(order_name, original_name))
                            # Update dict first
                            permanent_partial_orders[order_name] = order_id
                            if original_name in partial_orders:
                                permanent_partial_orders[original_name] = partial_orders[original_name]
                                partial_orders.pop(original_name, None)
                                rewrite_name_id_file(partial_orders_file, partial_orders)
                            # Atomically rewrite permanent_partial file
                            rewrite_name_id_file(permanent_partial_orders_file, permanent_partial_orders)
                        else:
                            # Normal partial: record for retry on next run
                            with open(partial_orders_file, "a") as f:
                                f.write("{},{}\n".format(order_name, order_id))
                            partial_orders[order_name] = order_id
                        
                        placed_orders.pop(order_name, None)
                        rewrite_name_id_file(placed_orders_file, placed_orders)
                        
                        if(failed_count != 0):
                            order_list_failed_download.append(num)
                        total_failed += failed_count
                        total_skipped += skipped_count
                        total_success += success_count
                        continue
                    
                    if state != "success":
                        print("\n Order not success, but status {}".format(state))
                        order_list_failed.append(num)
                        continue
                    
                    # SUCCESS
                    results = response['_links']['results']
                    
                    print("\n Total files to download: {} ".format(len(results)))
                    
                    total_files_to_download += len(results)
                    
                    filename = "order_result_{}.json".format(order_name)
                    with open(os.path.join(output_site_dir, filename), 'w') as outfile:
                        json.dump(response, outfile)
        
                    failed_count, skipped_count, success_count = download_results(results, output_site_dir, order_name)
                    
                    if(failed_count == 0):
                        order_id = placed_orders.get(order_name) or order_url.rsplit("/", 1)[-1]
                        with open(completed_orders_file, "a") as f:
                            f.write("{},{}\n".format(order_name, order_id))
                        completed_orders[order_name] = order_id
                        placed_orders.pop(order_name, None)
                        rewrite_name_id_file(placed_orders_file, placed_orders)
                        
                        # If this is a retry order, promote original partial to completed
                        if order_name.endswith("_retry"):
                            original_name = order_name[:-len("_retry")]
                            if original_name in partial_orders:
                                original_id = partial_orders[original_name]
                                with open(completed_orders_file, "a") as f:
                                    f.write("{},{}\n".format(original_name, original_id))
                                completed_orders[original_name] = original_id
                                partial_orders.pop(original_name, None)
                                rewrite_name_id_file(partial_orders_file, partial_orders)
                                print(" Retry succeeded - original {} promoted to completed".format(original_name))
                    
                    if(failed_count != 0):
                        order_list_failed_download.append(num)
                    
                    total_failed += failed_count
                    total_skipped += skipped_count
                    total_success += success_count
            
                except Exception as e:
                    print("Error occurred:", str(e))
                    pass
            
            print("\nSUMMARY")
            print("Total Files to be downloaded: {}".format(total_files_to_download))
            print("Total Files Failed to Download: {}".format(total_failed))
            print("Total Files Skipped: {}".format(total_skipped))
            print("Total Files Success: {}".format(total_success))
            print("Order Chunks Failed: {}".format(order_list_failed))
            print("Download Chunks Failed: {}".format(order_list_failed_download))
            print("--- %s seconds ---" % (time.time() - start_time))

print(sys.argv)

if __name__ == "__main__":
    if "--parallel" in sys.argv:
        multiprocessing.set_start_method('fork', force=True)
    main(sys.argv[1:])
