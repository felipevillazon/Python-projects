import pandas as pd
import sys
import os

# This function will create an ini file from a file.xlsx
# The input parameters are the path to the xlsx file and the location you would like to place the ini file


# This function compares if the headers of the xlsx file match the expected headers
def compare_headers(expected, actual):
    # Initialize a list to store missing headers
    missing_headers = []

    # Check if each expected header is in the actual headers
    for header in expected:
        if header not in actual:
            missing_headers.append(header)

    # If there are missing headers, return the list of missing headers
    if missing_headers:
        return False, missing_headers
    
    # If no headers are missing, return True and an empty list
    return True, missing_headers
    
    
# Function to check if all mandatory headers are present in the actual headers
def check_mandatory_headers(mandatory, actual):
    missing_headers = [header for header in mandatory if header not in actual]
    if missing_headers:
        print("The following mandatory headers are missing:")
        #for header in missing_headers:
         #   print(f"- {header}")
        return missing_headers
    return True


def create_ini_file_from_excel(excel_file, ini_file_path):
    
    # Expected headers in the xlsx file
    expected_headers = ['number', 'p&id number', 'p&id tag', 'sensor PID tag', 'Device','Equipment' , 'description', 'function description', 'Name', 'Data Type', 'Start value',
                                  'Retain', 'Accessible ', 'Writable', 'Visible', 'Setpoint', 'Comment', 'Log']
                                  
    mandatory_headers = ['p&id number', 'Name', 'Data Type', 'description', 'Log']
                        
    # Define the path for the error log file
    error_file_path = os.path.splitext(ini_file_path)[0] + '_errors.log'
    
    # Read all sheets from the xlsx file
    sheets = pd.read_excel(excel_file, sheet_name=None, dtype=str)
    
    with open(ini_file_path, 'w') as f, open(error_file_path, 'w') as error_file:
    
        
        # We loop over all the different sheets present in the xlsx file
        for sheet_name, df in sheets.items():
        
          
            actual_headers = df.columns.tolist()
            if actual_headers != expected_headers:
                missing_headers = compare_headers(expected_headers, actual_headers)
                missing_mandatory_headers = check_mandatory_headers(mandatory_headers,actual_headers)
                print(missing_mandatory_headers)
                error_file.write("\n")
                error_file.write(f"Processing sheet format...\n")
                error_file.write("\n")
                
                error_file.write(f"Warning: You are missing a column in sheet '{sheet_name}', the headers are: '{missing_headers[1]}'\n")
                if len(missing_mandatory_headers)!=0:
                    # Raise error when finding a difference in the headers. Program stops.
                    raise ValueError(f"Sheet '{sheet_name}' does not have the expected headers. You are missing a mandatory header: '{missing_mandatory_headers}'")
            
            # Write section header with sheet name
            f.write(f"[{sheet_name}]\n")
            # We leave a space between the header and the data
            f.write(f"\n")
            
            error_file.write("\n")
            error_file.write(f"Processing sheet '{sheet_name}'...\n")
            error_file.write("\n")
            # loop over the different rows
            for index, row in df.iterrows():
                    
                    # Looking for missing values. The error aare saved in a file as a warning message.
                    
                    if pd.isna(row['p&id number']):
                        error_file.write(f"Warning: NaN value found in 'p&id number' in sheet '{sheet_name}', row {index + 1}\n")
                        continue
                        
                    if pd.isna(row['Name']):
                        error_file.write(f"Warning: NaN value found in 'Name' in sheet '{sheet_name}', row {index + 1}\n")
                    
                    if pd.isna(row['Data Type']):
                        error_file.write(f"Warning: NaN value found in 'Data Type' in sheet '{sheet_name}', row {index + 1}\n")
                        
                    if pd.isna(row['description']):
                        error_file.write(f"Warning: NaN value found in 'description' in sheet '{sheet_name}', row {index + 1}\n")
                        
                    if pd.isna(row['Log']):
                        error_file.write(f"Warning: NaN value found in 'Log' in sheet '{sheet_name}', row {index + 1}\n")
                        
                    
                    # Write in the ini file the selected items from the xlsx file
                    try:
                        f.write(f"    db_num={str(row['p&id number']).zfill(3)}\n")
                        f.write(f"    [{sheet_name}.{row['Name']}]\n")
                        f.write(f"       data_type={row['Data Type']}\n")
                        f.write(f"       description={row['description']}\n")
                        f.write(f"    log={row['Log']}\n")
                        f.write(f"\n")
                    except Exception as e:
                        error_file.write(f"Error processing row {index + 1} in sheet '{sheet_name}': {e}\n")
            # Add an empty line after each section
            f.write("\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py input_file.xlsx output_file.ini")
        sys.exit(1)

    excel_file = sys.argv[1]
    ini_file_path = sys.argv[2]

    create_ini_file_from_excel(excel_file, ini_file_path)

    print(f"INI file '{ini_file_path}' has been created successfully.")

