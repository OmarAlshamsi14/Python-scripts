import tkinter as tk 
from tkinter import simpledialog
from tkinter import messagebox
import os
import csv

from Bio import Entrez # put the modules in file for this to work!!!!!!!!!!!!!!!!!!!!!!
from Bio import SeqIO

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


def sequenceGen():
    filename = "output_file.csv"
    f = open(filename, "r")
    header = f.readline().strip().split() #Header to ignore
    lines = f.readlines() #Get all lines (Besides header)
    f.close()

    o = open("Link_Output.csv", "a", newline='')
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
          
    log_text.insert(tk.END, "Done generating sequeneces\n")
    
    

def delete_csv_files():
    folder_path = os.path.dirname(os.path.realpath(__file__))
    csv_files = [file for file in os.listdir(folder_path) if 'rnaseq_data' in file]
    for file in csv_files:
        os.remove(os.path.join(folder_path, file))
    log_text.insert(tk.END, "CSV files have been deleted.\n")
    root.update_idletasks()
    
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
    product = (simpledialog.askstring("Product Parameter seperated by commas", "Enter the product parameter:")).split(",")
    mcheck = 0
    if product:
        for item in product:
            fileautomation(item)
        if var1.get() == 1:
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
