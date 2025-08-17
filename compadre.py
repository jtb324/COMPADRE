from xopen import xopen
import os, math, sys, json, signal, socket, traceback, random
import pandas as pd

import warnings
warnings.filterwarnings('ignore', category=FutureWarning, module='pandas')

# COMPADRE imports
import ersa
from pop_classifier import run_new as run_pop_classifier

########################################

def signal_handler(signum, frame):
    try:
        if 'server_socket' in globals():
            server_socket.close()
    except:
        pass
    os._exit(0)

def safe_print(*args, **kwargs):
    try:
        print(*args, **kwargs)
        sys.stdout.flush()
    except BrokenPipeError:
        # Python flushes standard streams on exit, so redirect remaining output to devnull to avoid another BrokenPipeError at shutdown
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(0)  # Python exits with error code 1 on EPIPE

def convert_value(value):

    """Convert string values of ERSA options to appropriate Python types"""
    # Try int first
    try:
        return int(value)
    except ValueError:
        pass
    
    # Try float
    try:
        return float(value)
    except ValueError:
        pass
    
    # Handle empty string as None for file paths
    if value == "None":
        return None
        
    # Return as string
    return value

def calculate_ersa_props(model_df):

    if model_df.empty:
        return 0, 0, 0, 1

    model_df['degree_of_relatedness'] = model_df['degree_of_relatedness'].astype(int)
    model_df['maxlnl2'] = model_df['maxlnl'].apply(lambda x: math.exp(x))
    
    model_df_2d = model_df[model_df['degree_of_relatedness'] == 2] 
    model_df_2d_sum = model_df_2d['maxlnl2'].sum()
    model_df_3d = model_df[model_df['degree_of_relatedness'] == 3] 
    model_df_3d_sum = model_df_3d['maxlnl2'].sum()

    ## Old version
    # model_df_4d = model_df[model_df['degree_of_relatedness'] == 4] 
    # model_df_4d_sum = model_df_4d['maxlnl2'].sum()
    # model_df_un = model_df[(model_df['degree_of_relatedness'] >= 5) & (model_df['degree_of_relatedness'] <= 40)] 
    # model_df_un_sum = model_df_un['maxlnl2'].sum()

    ## New version -- updated 7.2.24
    # model_df_4d = model_df[(model_df['degree_of_relatedness'] >= 4) & (model_df['degree_of_relatedness'] < 10)]
    # model_df_4d_sum = model_df_4d['maxlnl2'].sum()
    # model_df_un = model_df[(model_df['degree_of_relatedness'] >= 10) & (model_df['degree_of_relatedness'] <= 40)] 
    # model_df_un_sum = model_df_un['maxlnl2'].sum()

    # Newest version
    model_df_4d = model_df[(model_df['degree_of_relatedness'] >= 4) & (model_df['degree_of_relatedness'] < 8)]
    model_df_4d_sum = model_df_4d['maxlnl2'].sum()
    model_df_un = model_df[(model_df['degree_of_relatedness'] >= 8) & (model_df['degree_of_relatedness'] <= 40)] 
    model_df_un_sum = model_df_un['maxlnl2'].sum()

    # re-proportion them 
    total_prop = model_df_2d_sum + model_df_3d_sum + model_df_4d_sum + model_df_un_sum
    model_df_2d_prop = model_df_2d_sum/total_prop
    model_df_3d_prop = model_df_3d_sum/total_prop
    model_df_4d_prop = model_df_4d_sum/total_prop
    model_df_un_prop = model_df_un_sum/total_prop

    return model_df_2d_prop, model_df_3d_prop, model_df_4d_prop, model_df_un_prop

def clean_ersa_options(ersa_option_dict):

    # Check data types and revert to default values if they're not correct

    int_options = {
        'ascertained_position': -1,
        'mask_region_sim_count': 0,
        'mask_region_cross_length': 1000000,
        'mask_region_threshold': 4,
        'number_of_chromosomes': 22,
        'max_meioses': 40,
    }
    float_options = {
        'parent_offspring_zscore': 2.33,
        'pois_mean': 13.73,
        'exp_mean': 3.197036753,
        'rec_per_meioses': 35.2548101,
        'min_cm': 2.5,
        'confidence_level': 0.9,
        'max_cm': 10
    }

    return {k: v if isinstance(v, (int, float)) else (int_options.get(k, v) if k in int_options else float_options.get(k, v)) for k, v in ersa_option_dict.items()}


def find_available_port(preferred_port, host='localhost', max_attempts=100):
    
    # First try the preferred port
    test_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    test_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    try:
        test_socket.bind((host, preferred_port))
        test_socket.close()
        return preferred_port
    except OSError:
        test_socket.close()
        safe_print(f"Port {preferred_port} is in use, searching for available port...", file=sys.stderr)
    
    # If preferred port is taken, try random ports in a range
    for _ in range(max_attempts):
        port = random.randint(4001, 8000)  # Try ports in this range
        test_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        test_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            test_socket.bind((host, port))
            test_socket.close()
            return port
        except OSError:
            test_socket.close()
            continue
    
    raise Exception(f"Could not find an available port after {max_attempts} attempts")

########################################

def main(segment_data_file, portnumber, ersa_flag_str, output_directory):

    signal.signal(signal.SIGINT, signal_handler)   # CTRL+C
    signal.signal(signal.SIGTERM, signal_handler)  # Termination
    signal.signal(signal.SIGPIPE, signal_handler)  # Broken pipe
    global server_socket

    segment_dict = {}
    additional_options = {}
    ibd2_status = 'false'

    if segment_data_file != 'NA':
        
        # First step -- read in ERSA options from COMPADRE runtime flags
        ersa_flags = ersa_flag_str.strip('"').split('|')
        flag_names = ersa_flags[::2]
        flag_values = ersa_flags[1::2]
        
        for flag, value in zip(flag_names, flag_values):
            additional_options[flag] = convert_value(value)

        additional_options = clean_ersa_options(additional_options)  

        with xopen(segment_data_file, 'r') as f: # Open the large file here and populate dictionary that stays in system memory

            # Check number of columns to determine if it has IBD1/2 data or not 
            header_line = f.readline().strip()
            num_columns = len(header_line.split())

            if num_columns == 7: # File with segment-by-segment IBD status in 7th column 
                ibd2_status = 'true'

                if 'chrom' in header_line: #skip header if it exists
                    next(f) # skip header line in this file version

                for line in f:
                    ls = line.strip().split('\t')
                    iid1, iid2, start, end, cmlen, chrom, ibd = ls[0], ls[1], int(ls[2]), int(ls[3]), round(float(ls[4]), 2), int(ls[5]), int(ls[6].strip())
                    key = f"{iid1}:{iid2}"
                    value = (chrom, start, end, cmlen, ibd)

                    if cmlen >= 5.0: 
                        if key not in segment_dict:
                            segment_dict[key] = [value,]
                        else:
                            segment_dict[key] += [value,]

            elif num_columns == 6: # File without segment-by-segment IBD status 
                
                if 'chrom' in header_line: #skip header if it exists
                    next(f)

                for line in f:
                    ls = line.split('\t')
                    if len(ls) == 11: # germline1
                        iid1, iid2, start, end, cmlen, chrom = ls[0].split(' ')[0], ls[1].split(' ')[0], int(ls[3].split(' ')[0]), int(ls[3].split(' ')[1]), round(float(ls[6]), 2), int(ls[2])
                    else: # germline2
                        iid1, iid2, start, end, cmlen, chrom = ls[0], ls[1], int(ls[2]), int(ls[3]), round(float(ls[4]), 2), int(ls[5].strip())
                    key = f"{iid1}:{iid2}"
                    value = (chrom, start, end, cmlen, 'NA') # no ibd1/2 data in this input

                    if cmlen >= 5.0:
                        if key not in segment_dict:
                            segment_dict[key] = [value,]
                        else:
                            segment_dict[key] += [value,]

            else:
                safe_print('[COMPADRE] Unrecognized segment file format. Please refer to the README (https://github.com/belowlab/compadre) for formatting guidelines.')
                sys.exit(1)

    #########################

    socket_host = os.environ['COMPADRE_HOST'] if 'COMPADRE_HOST' in os.environ else 'localhost'

    # Find an available port
    actual_port = find_available_port(portnumber, socket_host)
    
    # If port changed, notify via stderr first (before stdout message)
    if actual_port != portnumber:
        safe_print(f"PORT_CHANGED:{actual_port}", file=sys.stderr)
    
    server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    server_socket.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)
    server_socket.settimeout(None)

    try:
        server_socket.bind((socket_host, actual_port))
        server_socket.listen(1)
        safe_print(f"COMPADRE helper is ready on port {actual_port}.")
        sys.stdout.flush()

        while True:
            try:
                conn, address = server_socket.accept()
                conn.settimeout(None)

                # NEW: Handle the client connection with better error handling
                handle_client_connection(conn, segment_dict, additional_options, ibd2_status, segment_data_file, output_directory)

            except KeyboardInterrupt:
                break
            except Exception as e:
                safe_print(f"Error accepting connection: {str(e)}")
                import traceback
                safe_print(f"Traceback: {traceback.format_exc()}")
                continue

    except KeyboardInterrupt:
        pass
    finally:
        try:
            server_socket.close()
        except:
            pass
        safe_print("[COMPADRE] COMPADRE helper shutdown complete.")


def handle_client_connection(conn, segment_dict, additional_options, ibd2_status, 
                           segment_data_file, output_directory):
    
    try:
        msg = conn.recv(1024).decode()
        
        if msg == 'close':
            conn.send("Closing server".encode())
            return
            
        if msg.strip() == 'test':
            conn.send("OK".encode())
            return

        ms = msg.strip().split('|')

        ########################################################
        # Population classifier
        if len(ms) >= 3 and ms[-1] == 'pop_classifier':
            try:
                #safe_print(f"Processing pop_classifier request", file=sys.stderr)
                
                eigenvec_file = ms[0]
                pop_file = ms[1]

                predictions = run_pop_classifier(eigenvec_file, pop_file)
                predictions = "|".join(str(x) for x in predictions)

                conn.send(predictions.encode())
                
            except Exception as e:
                error_msg = f"POP_CLASSIFIER_ERROR: {str(e)}\n{traceback.format_exc()}"
                safe_print(error_msg, file=sys.stderr)
                conn.send(f"ERROR: {error_msg}".encode())

        ########################################################
        # PADRE
        elif len(ms) >= 1 and ms[-1] == 'padre':
            try:
                safe_print(f"Processing PADRE request", file=sys.stderr)
                
                # object to pass to ersa is the WHOLE dictionary
                segment_obj = json.dumps(segment_dict)

                # make output directory
                ersa_dir = f'{output_directory}/ersa'
                if not os.path.exists(ersa_dir):
                    os.makedirs(ersa_dir, exist_ok=True)
                ersa_outfile = f'{ersa_dir}/output_all_ersa'
                
                ersa_options = {
                    "segment_dict": segment_obj,
                    "segment_files": segment_data_file,
                    "model_output_file": f"{ersa_outfile}.model",
                    "output_file": f"{ersa_outfile}.out",
                    "return_output": False,
                    "write_output": True,
                    "use_ibd2_siblings": ibd2_status 
                }

                # merge with additional options
                ersa_options.update(additional_options)

                #safe_print(f"Running ERSA with options: {ersa_options}", file=sys.stderr)
                ersa.runner(ersa_options) # output written to file, not returned

                conn.send(ersa_outfile.encode()) # return filepath it was written to
                
            except Exception as e:
                error_msg = f"PADRE_ERROR: {str(e)}\n{traceback.format_exc()}"
                safe_print(error_msg, file=sys.stderr)
                conn.send(f"ERROR: {error_msg}".encode())

        ########################################################
        # Pairwise ERSA processing
        elif len(ms) >= 4:
            try:
                id1, id2, vector_str, analysis_type = ms
                #safe_print(f"Processing ERSA for {id1} <-> {id2}", file=sys.stderr)

                # Process the pairwise request
                result = process_pairwise_ersa(id1, id2, vector_str, analysis_type, 
                                             segment_dict, additional_options, 
                                             ibd2_status, segment_data_file, output_directory)
                
                conn.send(result.encode())
                
            except Exception as e:
                error_msg = f"ERSA_ERROR: Failed to process {id1} <-> {id2}: {str(e)}\n{traceback.format_exc()}"
                safe_print(error_msg, file=sys.stderr)
                conn.send(f"ERROR: {error_msg}".encode())

        else:
            error_msg = f"UNKNOWN_COMMAND: Received message with {len(ms)} parts: {ms}"
            safe_print(error_msg, file=sys.stderr)
            conn.send(f"ERROR: {error_msg}".encode())

    except Exception as e:
        error_msg = f"SOCKET_ERROR: {str(e)}\n{traceback.format_exc()}"
        safe_print(error_msg, file=sys.stderr)
        try:
            conn.send(f"ERROR: {error_msg}".encode())
        except:
            pass  # Socket might be closed
    finally:
        try:
            conn.close()
        except:
            pass


def process_pairwise_ersa(id1, id2, vector_str, analysis_type, segment_dict, 
                         additional_options, ibd2_status, segment_data_file, output_directory):
    
    try:
        # Step 1: Parse inputs
        try:
            id1_temp = id1.split('_')[-1] # just in case
            id2_temp = id2.split('_')[-1]
            
            idcombo = f"{id1_temp}:{id2_temp}"
            idcombo2 = f"{id2_temp}:{id1_temp}"
            
            vector_arr = [float(x) for x in vector_str.split(',')]
            if len(vector_arr) != 6:
                raise ValueError(f"Expected 6 values in vector, got {len(vector_arr)}")
                
        except Exception as e:
            raise Exception(f"Failed to parse inputs (id1={id1}, id2={id2}, vector='{vector_str}'): {str(e)}")

        # Step 2: Check for segments
        if idcombo not in segment_dict and idcombo2 not in segment_dict:
            #safe_print(f"No segments found for {idcombo} or {idcombo2}", file=sys.stderr)
            return '0,0,0,0,0,1'
        
        # Step 3: Prepare segment data
        try:
            if idcombo in segment_dict:
                key = idcombo
                segment_obj = {key : segment_dict[key]} 
            else:
                key = idcombo2
                segment_obj = {key : segment_dict[key]} 
                
            segment_obj = json.dumps(segment_obj)
            #safe_print(f"Found {len(segment_dict[key])} segments for {key}", file=sys.stderr)
            
        except Exception as e:
            raise Exception(f"Failed to prepare segment data for {key}: {str(e)}")

        # Step 4: Set up ERSA options
        try:
            ersa_dir = f'{output_directory}/ersa'
            if not os.path.exists(ersa_dir):
                os.makedirs(ersa_dir, exist_ok=True)
            ersa_outfile = f'{ersa_dir}/output_all_ersa'

            ersa_options = {
                "single_pair": key,
                "segment_dict": segment_obj,
                "segment_files": segment_data_file,
                "model_output_file": f"{ersa_outfile}.model",
                "output_file": f"{ersa_outfile}.out",
                "return_output": True,
                "write_output": False,
                "use_ibd2_siblings": ibd2_status,
            }

            # merge with additional options
            ersa_options.update(additional_options)
            
        except Exception as e:
            raise Exception(f"Failed to set up ERSA options: {str(e)}")

        # Step 5: Run ERSA
        try:
            #safe_print(f"Running ERSA for {key}", file=sys.stderr)
            output_model_df = ersa.runner(ersa_options)
            
        # except Exception as e:
        #     safe_print(f"ERSA failed for {key} with segments: {segment_dict[key]}", file=sys.stderr)
        #     raise Exception(f"ERSA runner failed: {str(e)}")

        except Exception as ersa_error:
            # Log the ERSA error but don't crash - return original vector
            safe_print(f"ERSA failed for {key} with segments: {segment_dict[key]}, reverting to original vector: {str(ersa_error)}", file=sys.stderr)
            return vector_str

        # Step 6: Process results
        try:
            if output_model_df.empty:
                safe_print(f"ERSA returned empty dataframe for {key}", file=sys.stderr)
                return '0,0,0,0,0,1'

            if output_model_df['maxlnl'].isna().all():
                safe_print(f"ERSA returned all NaN maxlnl values for {key}", file=sys.stderr)
                return '0,0,0,0,0,1'
            
            ersa_props = calculate_ersa_props(output_model_df)
            prop02 = 1 - (vector_arr[0] + vector_arr[1])
            ersa_props_updated = tuple(x * prop02 for x in ersa_props) 
            updated_vector = f'{vector_arr[0]},{vector_arr[1]},{ersa_props_updated[0]},{ersa_props_updated[1]},{ersa_props_updated[2]},{ersa_props_updated[3]}'
            
            #safe_print(f"ERSA completed successfully for {key}", file=sys.stderr)
            return updated_vector
            
        except Exception as e:
            raise Exception(f"Failed to process ERSA results: {str(e)}")

    except Exception as e:
        # Re-raise with more context
        raise Exception(f"ERSA processing failed for {id1} <-> {id2}: {str(e)}")


if __name__ == '__main__':

    # This is how the script is called from bin/compadre_kickoff.pl

    segment_data_file = sys.argv[1]
    portnumber = int(sys.argv[2])
    ersa_flag_str = sys.argv[3]
    output_directory = sys.argv[4]

    main(segment_data_file, portnumber, ersa_flag_str, output_directory)