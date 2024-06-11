June 10, 2024
- Functions introduced: 
  1) trigger_enum: Provides a direct way to enumerate files from the UI.
  2) trigger_generate_fasta_file: Allows users to trigger the FASTA file generation directly from the UI.
  3) copy_and_rename_files: Manages copying and renaming output files based on product parameters, ensuring organized and clear results. Copying the file is a way to make the FASTA function actually work (still in progress) 

Significant effort has been made to improve the generate_fasta_file function. A helper method now copies the content of output_file.csv before renaming it, which is a key step in resolving the issue; however, the function is still under progress.
