import tkinter as tk
from tkinter import simpledialog, messagebox
import os
import shutil #for copying file
import csv
import concurrent.futures # for asynchronously executing
from Bio import Entrez
from Bio import SeqIO

# Set your email for Entrez
Entrez.email = "omaralshamsi14@gmail.com"  # Change this to your actual email address

def enum(): #Gives every rnaseq data file a number to later start filtering.
    folder_path = os.path.dirname(os.path.realpath(__file__))
    log_text.insert(tk.END, "Enumerating files...\n")
    root.update_idletasks()

    # Get a list of all CSV files in the folder
    csv_files = [file for file in os.listdir(folder_path) if file.startswith('rnaseq_data')]

    # Rename each file with a numbered format
    for index, file in enumerate(csv_files, start=1):
        old_path = os.path.join(folder_path, file)
        new_name = f"{index}.csv"
        new_path = os.path.join(folder_path, new_name)
        os.rename(old_path, new_path)

    log_text.insert(tk.END, "Files have been successfully renamed.\n")
    root.update_idletasks()

def filter_and_append_to_csv(input_file, output_file, product):
    with open(input_file, "r") as input_csv, open(output_file, "a", newline='') as output_csv:
        csv_reader = csv.reader(input_csv)
        csv_writer = csv.writer(output_csv)

        # If the output file is empty, write the header to the output CSV file
        if output_csv.tell() == 0:
            header = next(csv_reader, None)
            if header:
                csv_writer.writerow(header)

        # Iterate through each row in the input CSV file
        for row in csv_reader:
            if len(row) > 0 and product in row[0]:
                # Write the matching row to the output CSV file
                csv_writer.writerow(row)

    log_text.insert(tk.END, f"Filtered and appended data to {output_file}.\n")
    root.update_idletasks()

def linkgen():
    filename = "output_file.csv"
    f = open(filename, "r")
    header = f.readline().strip().split() #Header to ignore
    lines = f.readlines() #Get all lines (Besides header)
    f.close()

    with open("Link_Output.csv", "a", newline='') as o:
        csv_writer = csv.writer(o)
        csv_writer.writerow(["Gene seq length", "Link to JGI page"])
        for line in lines:
            line = line.strip().split()
            link = "https://img.jgi.doe.gov/cgi-bin/mer/main.cgi?section=MetaGeneDetail&page=genePageMainFaa&taxon_oid=" + line[0] + "&data_type=assembled&gene_oid=" + line[2].split(",")[0]
            try:
                if type(int(line[9])) != int:
                    continue
                else:
                    mlist = [line[9], link]
                csv_writer.writerow(mlist)
            except:
                continue
    log_text.insert(tk.END, "Done generating links.\n")
    root.update_idletasks()

# The function uses the Entrez module to fetch amino acid sequences from the protein database
def fetch_amino_acid_sequence(gene_id):
    try:
        handle = Entrez.efetch(db="protein", id=gene_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return gene_id, str(record.seq)
    except Exception as e:
        return gene_id, f"Error: {e}"

# gets sequences and write them to a fasta file
def generate_fasta_file():
    csv_filename = "output_file.csv"
    fasta_filename = "amino_acid_sequences.fasta"

    log_text.delete(1.0, tk.END)

    if not os.path.exists(csv_filename):
        log_text.insert(tk.END, f"Error: {csv_filename} does not exist.\n")
        root.update_idletasks()
        return

    log_text.insert(tk.END, f"Reading gene IDs from {csv_filename}...\n")
    root.update_idletasks()

    gene_ids = []
    with open(csv_filename, "r") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        header = next(csv_reader)

        log_text.insert(tk.END, f"CSV Header: {header}\n")

        for row in csv_reader:
            log_text.insert(tk.END, f"Processing row: {row}\n")
            if len(row) >= 4:  # Assuming the gene ID is in the fourth column (index 3)
                gene_id = row[3]  # Modify this index based on the actual structure
                if gene_id.strip():  # Check if gene_id is not empty or just whitespace
                    gene_ids.append(gene_id)
                else:
                    log_text.insert(tk.END, f"Skipping row with empty gene ID: {row}\n")
            else:
                log_text.insert(tk.END, f"Skipping row with insufficient columns: {row}\n")

    if not gene_ids:
        log_text.insert(tk.END, "No gene IDs found. Ensure the input file has the correct data.\n")
        root.update_idletasks()
        return

    log_text.insert(tk.END, f"Found {len(gene_ids)} gene IDs.\n")
    root.update_idletasks()

    log_text.insert(tk.END, f"Generating FASTA file: {fasta_filename}\n")
    root.update_idletasks()

    with open(fasta_filename, "w") as fasta_file:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            future_to_gene_id = {executor.submit(fetch_amino_acid_sequence, gene_id): gene_id for gene_id in gene_ids}
            
            for future in concurrent.futures.as_completed(future_to_gene_id):
                gene_id = future_to_gene_id[future]
                try:
                    gene_id, sequence = future.result()
                    if "Error:" in sequence:
                        log_text.insert(tk.END, f"Error fetching sequence for {gene_id}: {sequence}\n")
                    else:
                        if sequence:
                            fasta_file.write(f">{gene_id}\n{sequence}\n")
                            log_text.insert(tk.END, f"Successfully fetched sequence for {gene_id}\n")
                        else:
                            log_text.insert(tk.END, f"No sequence returned for {gene_id}\n")
                except Exception as e:
                    log_text.insert(tk.END, f"Error processing sequence for {gene_id}: {e}\n")
                root.update_idletasks()

    log_text.insert(tk.END, "FASTA file generation complete.\n")
    root.update_idletasks()


def delete_csv_files():
    """Delete all enumerated CSV files."""
    folder_path = os.path.dirname(os.path.realpath(__file__))
    csv_files = [file for file in os.listdir(folder_path) if has_numbers(file)]
    for file in csv_files:
        os.remove(os.path.join(folder_path, file))
        
    log_text.insert(tk.END, "CSV files have been deleted.\n")
    root.update_idletasks()

def has_numbers(inputString):
    return any(char.isdigit() for char in inputString)

def copy_and_rename_files(product):
    folder_path = os.path.dirname(os.path.realpath(__file__))
    original_file = os.path.join(folder_path, "output_file.csv")
    product_file = os.path.join(folder_path, f"{product}_output.csv")

    if os.path.exists(original_file):
        shutil.copyfile(original_file, product_file)
        log_text.insert(tk.END, f"Copied and renamed to {product_file}.\n")
    else:
        log_text.insert(tk.END, f"Error: {original_file} does not exist.\n")
    
    # Renaming Link_Output.csv based on product parameter
    link_output_file = os.path.join(folder_path, "Link_Output.csv")
    renamed_link_output_file = os.path.join(folder_path, f"{product}_Link_Output.csv")

    if os.path.exists(link_output_file):
        shutil.copyfile(link_output_file, renamed_link_output_file)
        log_text.insert(tk.END, f"Copied and renamed to {renamed_link_output_file}.\n")
    else:
        log_text.insert(tk.END, f"Error: {link_output_file} does not exist.\n")

    root.update_idletasks()

def fileautomation(product):
    log_text.delete(1.0, tk.END)
    try:
        enum()
        try:
            i = 1
            while True:
                input_file_path = f"{i}.csv"
                if os.path.exists(input_file_path):
                    filter_and_append_to_csv(input_file_path, "output_file.csv", product)
                else:
                    break
                i += 1
        except Exception as e:
            log_text.insert(tk.END, f"Error in filtering process: {e}\n")
            root.update_idletasks()
        
        # Generate links and FASTA file before renaming the output file
        linkgen()
        generate_fasta_file()

        # Copy and rename the output files after processing
        copy_and_rename_files(product)
        log_text.insert(tk.END, "Process completed successfully!\n")
    except Exception as e:
        messagebox.showinfo("Error", f"Error: {e}\n")
        log_text.insert(tk.END, f"Error: {e}\n")
    root.update_idletasks()

def get_product_parameter(): #For UI and getting product name
    root = tk.Tk()
    root.withdraw()  # Hide the main window
    product = simpledialog.askstring("Product Parameter", "Enter the product parameter:")
    if product:
        fileautomation(product)
    
    else:
        messagebox.showwarning("File Automation", "Product parameter not provided!")

def get_multi_product_parameter(): #For UI and getting product name
    root = tk.Tk()
    root.withdraw()  # Hide the main window
    product = (simpledialog.askstring("Product Parameter seperated by commas", "Enter the product parameter:"))
    if product is None:
         messagebox.showwarning("File Automation", "Product parameter not provided!")
    else:
        product = product.split(",")
    
    mcheck = 0
    if product:
        for item in product:
            fileautomation(item)
       
        delete_csv_files()
    else:
        messagebox.showwarning("File Automation", "Product parameter not provided!")


# Functions to directly trigger enum and generate_fasta_file when we hit the button
def trigger_enum():
    """Trigger the enumeration of files."""
    enum()

def trigger_generate_fasta_file():
    """Trigger the generation of the FASTA file."""
    generate_fasta_file()

if __name__ == "__main__":
    root = tk.Tk()
    root.title("File Automation")

    frame = tk.Frame(root)
    frame.pack(padx=10, pady=10)

    log_text = tk.Text(frame, height=10, width=50)
    log_text.pack(side=tk.TOP, padx=5, pady=5)

    button1 = tk.Button(frame, text="Enter Single Product Parameter", command=get_product_parameter)
    button1.pack(side=tk.TOP, padx=5, pady=5)

    button2 = tk.Button(frame, text="Enter Multiple Product Parameters", command=get_multi_product_parameter)
    button2.pack(side=tk.TOP, padx=5, pady=5)

    button3 = tk.Button(frame, text="Delete CSV Files", command=delete_csv_files)
    button3.pack(side=tk.TOP, padx=5, pady=5)

    button4 = tk.Button(frame, text="Enumerate Files", command=trigger_enum)
    button4.pack(side=tk.TOP, padx=5, pady=5)

    button5 = tk.Button(frame, text="Generate FASTA File", command=trigger_generate_fasta_file)
    button5.pack(side=tk.TOP, padx=5, pady=5)

    root.mainloop()
