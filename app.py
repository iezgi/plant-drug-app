from flask import Flask, render_template, request
from drug_utils import find_plants_for_disease

app = Flask(__name__)

# File paths
DRUG_FILE = "drug_data_example_FINAL.csv"
PLANT_FILE = "plant_chemicals_smiles_example_FINAL.csv"


@app.route("/", methods=["GET", "POST"])
def index():
    results = None
    if request.method == "POST":
        disease_code = request.form.get("disease")
        region = request.form.get("region")

        # Run your function
        results = find_plants_for_disease(DRUG_FILE, PLANT_FILE, disease_code, region)

    return render_template("index.html", results=results)


if __name__ == "__main__":
    app.run(debug=True)

