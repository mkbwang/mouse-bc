
import requests
import xml.etree.ElementTree as ET
import json
import re



def ncbi_search(ID):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    search_url = f"{base_url}esearch.fcgi?db=nuccore&term={ID}&retmode=xml"

    # Make the request to the NCBI API
    search_response = requests.get(search_url)

    # Check if the request was successful
    if search_response.status_code == 200:
        # Parse the XML response to find the UID
        root = ET.fromstring(search_response.content)
        uid_element = root.find(".//Id")
        if uid_element is not None:
            uid = uid_element.text
            # --- Step 2: Use the UID with esummary to get the gene name ---
            summary_url = f"{base_url}esummary.fcgi?db=nuccore&id={uid}&retmode=json"
            
            summary_response = requests.get(summary_url)
            
            if summary_response.status_code == 200:
                # Parse the JSON response
                summary_data = summary_response.json()
                # Navigate the JSON structure to find the gene symbol
                gene_info = summary_data['result'][uid]['title'].split(',')[0]
                return(gene_info)
            else:
                return f"{summary_response.status_code}"
        else:
            return "No gene symbol"
    else:
        return f"{search_response.status_code}"



def ensembl_search(ID):

    # Define the Ensembl REST API endpoint
    server = "https://rest.ensembl.org"
    endpoint = f"/xrefs/id/{ID}?"

    # Construct the full URL
    url = server + endpoint

    try:
        # Make the GET request
        response = requests.get(url, headers={"Content-Type": "application/json"})
        response.raise_for_status()  # Raise an error for bad status codes

        # Decode the JSON response
        data = response.json()

        # Find the gene symbol in the response
        gene_symbol = None
        for item in data:
            if item.get("info_type") == "MISC":
                gene_symbol = item.get("display_id")
                break

        if gene_symbol:
            gene_symbol = gene_symbol.split("-")[0]
            return gene_symbol
        else:
            return "No gene symbol"

    except:
        return "Error Fetching"

