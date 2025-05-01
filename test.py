from pybiomart import Server

server = Server(host="http://www.ensembl.org")
mart = server.marts['ENSEMBL_MART_SNP']  # Connect to the Variation mart
dataset = mart.datasets['hsapiens_snp']  # Access the SNP dataset

# Query SNP data
snp_data = dataset.query(
    attributes=['refsnp_id', 'chromosome_name', 'start_position', 'end_position'],
    filters={'refsnp_id': ['rs123456']}
)
print(snp_data)