import pandas as pd
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

class DrugData:
    def __init__(self, filepath):
        self.df = pd.read_csv(filepath)
        self.df.rename(columns={
            'Drug Name': 'Drug_Name',
            'ICD-10 Code': 'ICD10_Code'
        }, inplace=True)

    def get_chemical_by_disease(self, disease_name):
        matches = self.df[self.df['ICD10_Code'].str.contains(disease_name, case=False, na=False)]
        return matches[['Drug_Name', 'SMILES']].drop_duplicates()


class PlantData:
    def __init__(self, filepath):
        self.df = pd.read_csv(filepath)
        self.df.rename(columns={
            'Plant Name': 'Plant',
            'Country of Origin': 'Country_of_Origin'
        }, inplace=True)

        self.fingerprint_generator = GetMorganGenerator(radius=2, fpSize=2048)

    def get_plants_by_region(self, region_name):
        return self.df[self.df['Country_of_Origin'].str.contains(region_name, case=False, na=False)]

    def smiles_to_mol(self, smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"[WARNING] Could not parse SMILES: {smiles}")
            return mol
        except Exception as e:
            print(f"[ERROR] Invalid SMILES '{smiles}': {e}")
            return None

    def compute_fingerprint(self, mol):
        return self.fingerprint_generator.GetFingerprint(mol)

    def find_matching_plants(self, drug_chemicals, region_name, similarity_threshold=0.3):
        plants_in_region = self.get_plants_by_region(region_name)
        matching_plants = []

        drug_fps = []
        for _, drug_row in drug_chemicals.iterrows():
            drug_mol = self.smiles_to_mol(drug_row['SMILES'])
            if drug_mol:
                drug_fps.append(self.compute_fingerprint(drug_mol))

        for _, plant_row in plants_in_region.iterrows():
            plant_mol = self.smiles_to_mol(plant_row['SMILES'])
            if not plant_mol:
                continue

            plant_fp = self.compute_fingerprint(plant_mol)

            for drug_fp in drug_fps:
                similarity = DataStructs.TanimotoSimilarity(plant_fp, drug_fp)
                if similarity >= similarity_threshold:
                    matching_plants.append(plant_row['Plant'])
                    break

        return list(set(matching_plants))


def find_plants_for_disease(drug_file, plant_file, disease_name, region_name):
    drug_data = DrugData(drug_file)
    plant_data = PlantData(plant_file)

    drug_chemicals = drug_data.get_chemical_by_disease(disease_name)
    matching_plants = plant_data.find_matching_plants(drug_chemicals, region_name)

    return matching_plants
