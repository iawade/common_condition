import argparse
from bioservices import BioMart
import requests
import sys
import time
import logging

# Configure logging to capture warnings from bioservices
logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger('bioservices')

def get_gene_coordinates(ensembl_id, assembly='GRCh38',  max_retries=5, retry_delay=5):
    server = BioMart()

    # Set the dataset and attributes
    dataset = 'hsapiens_gene_ensembl'
    server.add_dataset_to_xml(dataset)

    # Add attributes for gene start and end positions
    server.add_attribute_to_xml("ensembl_gene_id")
    server.add_attribute_to_xml("chromosome_name")
    server.add_attribute_to_xml("start_position")
    server.add_attribute_to_xml("end_position")

    # Add filter for the Ensembl gene ID
    server.add_filter_to_xml("ensembl_gene_id", ensembl_id)

# Retry logic for handling Internal Server Errors or other exceptions
    attempt = 0
    while attempt < max_retries:
        try:
            # Execute the query
            xml_query = server.get_xml()
            result = server.query(xmlq=xml_query)

# If result contains warning or error, log and retry
            if isinstance(result, str) and "Internal Server Error" in result:
                logger.warning(f"Warning from BioMart: {result}")
                attempt += 1
                if attempt < max_retries:
                    print(f"Retrying... ({attempt}/{max_retries})")
                    time.sleep(retry_delay)  # Wait before retrying
                else:
                    print("Max retries reached. Exiting.", file=sys.stderr)
                    sys.exit(1)

            # If we get a valid result, break the loop
            elif result:
                break
            else:
                raise ValueError("Received empty result from BioMart")

        except requests.exceptions.RequestException as e:
            # Catch request exceptions like connection errors
            print(f"Error: Unable to connect to BioMart. The server may be down. Details: {e}", file=sys.stderr)
            attempt += 1
            if attempt < max_retries:
                print(f"Retrying... ({attempt}/{max_retries})")
                time.sleep(retry_delay)  # Wait before retrying
            else:
                print("Max retries reached. Exiting.", file=sys.stderr)
                sys.exit(1)
        except Exception as e:
            # Catch other exceptions
            print(f"An unexpected error occurred: {e}", file=sys.stderr)
            sys.exit(1)

    # Dynamically create the output filename
    output_file = f"{ensembl_id}.bed"

    # Process the result and write to the output file in BED format
    with open(output_file, 'w') as output:
        for line in result.split('\n'):
            if line:
                row = line.split('\t')
                gene_id = row[0]
                chromosome = row[1]
                start_position = row[2]
                end_position = row[3]

                output.write(f"chr{chromosome}\t{start_position}\t{end_position}\t{gene_id}\n")

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Retrieve gene coordinates from BioMart.')
    parser.add_argument('--ensembl_id', required=True, help='Ensembl gene ID')
    args = parser.parse_args()

    # Run the function with provided arguments
    get_gene_coordinates(args.ensembl_id)