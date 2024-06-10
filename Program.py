import tkinter as tk 
from tkinter import simpledialog
from tkinter import messagebox
import os
import csv
import concurrent.futures # for asynchronously executing
from Bio import Entrez 
from Bio import SeqIO

Entrez.email = "YourEmail@example.com"  # Change it to your email address

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
            header = next(csv_reader)
            csv_writer.writerow(header)

        # Iterate through each row in the input CSV file
        for row in csv_reader:
            if product in row[0]:
                # Write the matching row to the output CSV file
                csv_writer.writerow(row)

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
    handle = Entrez.efetch(db="protein", id=gene_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")  # Reading the sequence
    handle.close()
    return gene_id, str(record.seq)  # Returns a tuple with the gene_id and the amino acid sequence

# gets sequences and write them to a fasta file
def generate_fasta_file():
    filename = "output_file.csv"
    if not os.path.exists(filename):
        log_text.insert(tk.END, f"Error: {filename} does not exist.\n")
        root.update_idletasks()
        return

    with open(filename, "r") as f:
        header = f.readline().strip().split()
        lines = f.readlines()

    # Extracting the gene ID
    gene_ids = [line.strip().split()[2].split(",")[0] for line in lines]

    with open("amino_acid_sequences.fasta", "w") as fasta_file:
        # Uses ThreadPool for parallel processing
        with concurrent.futures.ThreadPoolExecutor() as executor:
            future_to_gene_id = {executor.submit(fetch_amino_acid_sequence, gene_id): gene_id for gene_id in gene_ids}
            for future in concurrent.futures.as_completed(future_to_gene_id):
                gene_id = future_to_gene_id[future]
                try:
                    _, sequence = future.result()
                    fasta_file.write(f">{gene_id}\n{sequence}\n")
                except Exception as e:
                    log_text.insert(tk.END, f"Error fetching sequence for {gene_id}: {e}\n")
                    root.update_idletasks()

    log_text.insert(tk.END, "Done generating FASTA file.\n")
    root.update_idletasks()

def sequenceGen():
    filename = "output_file.csv"
    f = open(filename, "r")
    header = f.readline().strip().split() #Header to ignore
    lines = f.readlines() #Get all lines (Besides header)
    f.close()

    o = open("Sequence_Output.csv", "a", newline='')
    csv_writer = csv.writer(o)
    csv_writer.writerow(["Gene seq length", "Amino Acid Sequence"])
    for line in lines:
        line = line.strip().split()
        aminoSequence = line[0] + line[2].split(",")[0] + "yayforNow"
        try:
            if type(int(line[9])) != int:
               continue
            else:
                mlist = [line[9], aminoSequence]
                csv_writer.writerow(mlist)
        except:
            continue
    
            # taxon id = line[0] +geneid =  line[2].split(",")[0]
            #record.format("fasta")
          
    log_text.insert(tk.END, "Done generating sequeneces\n")

def delete_csv_files():
    folder_path = os.path.dirname(os.path.realpath(__file__))
    csv_files = [file for file in os.listdir(folder_path) if has_numbers(file)]
    for file in csv_files:
        os.remove(os.path.join(folder_path, file))
        
    log_text.insert(tk.END, "CSV files have been deleted.\n")
    root.update_idletasks()


def has_numbers(inputString):
    return any(char.isdigit() for char in inputString)
    
def renamer(product): #Updates names (For now doesn't change anything)
    folder_path = os.path.dirname(os.path.realpath(__file__))

    # Get a list of all CSV files in the folder
    csv_files = [file for file in os.listdir(folder_path) if file.endswith('.csv')]

    # Rename each file with a numbered format
    for index, file in enumerate(csv_files, start=1):
        old_path = os.path.join(folder_path, file)
        if file == "output_file.csv":
            new_name = product + "output.csv" #Change name to new name or leave as default
            new_path = os.path.join(folder_path, new_name)
            os.rename(old_path, new_path)
        elif file== "Link_Output.csv":
            new_name = product + "Link_Output.csv" #Change name to new nameor leave as default
            new_path = os.path.join(folder_path, new_name)
            os.rename(old_path, new_path)

    log_text.insert(tk.END, "Output files have been successfully renamed.\n")
    root.update_idletasks()

def fileautomation(product):
    mcheck = 1
    log_text.delete(1.0, tk.END)  # Clear previous log
    try:
        enum()
        try:
            i = 1;
            while(True):
                input_file_path = str(i) + ".csv"
                filter_and_append_to_csv(input_file_path, "output_file.csv", product)
                i = i + 1;
        except:
            pass
        # if(var1.get() == 1 & mcheck == 1):
            #delete_csv_files() #Delets all csv files
        linkgen() #Generates link to JGI pages for each data set with wanted product
        renamer(product)
        log_text.insert(tk.END, "Process completed successfully!\n")
    except Exception as e:
        messagebox.showinfo(tk.END, f"Error: {e}\n")
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

if __name__ == "__main__": #UI
    root = tk.Tk()
    root.title("File Automation")
    
    frame = tk.Frame(root)
    frame.pack(padx=10, pady=10)

    log_text = tk.Text(frame, height=10, width=50)
    log_text.pack(side=tk.TOP, padx=5, pady=5)

    button1 = tk.Button(frame, text="Enter Single Product Parameter", command=get_product_parameter)
    button1.pack(side=tk.BOTTOM, padx=5, pady=5)

    button2 = tk.Button(frame, text="Enter Multiple Product Parameters", command=get_multi_product_parameter)
    button2.pack(side=tk.BOTTOM, padx=5, pady=5)

    button3 = tk.Button(frame, text="delete csv files", command=delete_csv_files)
    button3.pack(side=tk.BOTTOM, padx=5, pady=5)

    root.mainloop()
