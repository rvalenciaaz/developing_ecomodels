import sqlite3
import pandas as pd
import sys

def add_bacteria(cursor, species_id, genus, species, strain, taxonomy, genome_filename):
    # Insert a new bacterial species into the table
    cursor.execute("""
    INSERT INTO bacterial_species (species_id, genus, species, strain, taxonomy, genome_filename)
    VALUES (?, ?, ?, ?, ?, ?);
    """, (species_id, genus, species, strain, taxonomy, genome_filename))

if __name__ == "__main__":
    # Check if the CSV filename was passed as an argument
    if len(sys.argv) != 2:
        print("Usage: python add_bacteria_from_table.py <csv_filename>")
        sys.exit(1)
    
    # Get the CSV filename from the command-line arguments
    csv_filename = sys.argv[1]

    # Connect to SQLite database
    conn = sqlite3.connect('simulations.db')
    cursor = conn.cursor()

    # Read bacteria data from the CSV file using pandas
    try:
        df = pd.read_csv(csv_filename)
    except FileNotFoundError:
        print(f"Error: File '{csv_filename}' not found.")
        sys.exit(1)

    # Ensure that the dataframe has the necessary columns including 'id'
    required_columns = ['id', 'genus', 'species', 'strain', 'taxonomy', 'genome_filename']
    if not all(col in df.columns for col in required_columns):
        raise ValueError(f"CSV file is missing one or more required columns: {required_columns}")

    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        species_id = row['id']  # Use the 'id' column from the CSV as the species_id
        genus, species, strain, taxonomy, genome_filename = row['genus'], row['species'], row['strain'], row['taxonomy'], row['genome_filename']

        # Add bacteria to the database
        add_bacteria(cursor, species_id, genus, species, strain, taxonomy, genome_filename)

    # Commit changes and close connection
    conn.commit()
    cursor.close()
    conn.close()