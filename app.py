from flask import Flask, render_template, url_for, request, redirect
from mordred._base.descriptor import Descriptor
from similarity_calc import DataLoader, Similarity
import predictor
from rdkit import Chem
from rdkit.Chem import Draw
from thetaplots import calculate_values, plotting
import io
import random
from flask import Response
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import webbrowser

from base64 import b64encode
from io import BytesIO
import os

# tmpl_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'templates')
# stat_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static')

# app = Flask(__name__, template_folder=tmpl_dir, static_folder=stat_dir)
app = Flask(__name__)

@app.after_request
def add_header(r):
    """
    Add headers to both force latest IE rendering engine or Chrome Frame,
    and also to cache the rendered page for 10 minutes.
    """
    r.headers["Cache-Control"] = "no-cache, no-store, must-revalidate"
    r.headers["Pragma"] = "no-cache"
    r.headers["Expires"] = "0"
    r.headers['Cache-Control'] = 'public, max-age=0'
    return r


@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        return render_template('similarity.html')
    return render_template('index.html')

@app.route('/similarity', methods=['POST', 'GET'])
def similarity():
    if request.method == 'POST':
        user_smiles = request.form['content']
        comparing_dataset = DataLoader.csv_loader('data/osn_data_similarity.csv')
        try:
            similar = Similarity.rejection_similarity(comparing_dataset, reference=user_smiles)
        except:
            return 'Wrong SMILES, please try again'

        Draw.MolToFile(similar[1],'static/dataset_image.png') 
        retrieved_features = [user_smiles, 
                                str(similar[2].iloc[0]),
                                round(similar[0][1],2), 
                                int(similar[2].iloc[1]*100),
                                int(similar[2].iloc[2]*100),
                                int(similar[2].iloc[3]*100),
                                int(similar[2].iloc[4]*100),
                                int(similar[2].iloc[5]*100),
                                int(similar[2].iloc[6]*100)]
        try:
            return render_template("similarity.html", tasks = retrieved_features)
        except:
            return 'There was an error adding your molecule. Please check the SMILES if they are correct. The model works reliably for molecules between 50-1000 g/mol having no metal ions.'

    if request.method == "GET":
        return render_template("similarity.html")


@app.route('/pls-prediction', methods=['POST', 'GET'])
def pls_prediction():
    if request.method == 'POST':
        user_smiles = request.form['content']
        #comparing_dataset = DataLoader.csv_loader('data/osn_data_similarity.csv')
        try:
        #    similar = Similarity.rejection_similarity(comparing_dataset, reference=user_smiles)
            user_descr = predictor.descripter(user_smiles)
            prediction = predictor.predictor(user_descr)
        except:
            return 'Wrong SMILES, please try again. Currently, molecules with more than 5 heavy atoms work only.'

        Draw.MolToFile(Chem.MolFromSmiles(user_smiles),'static/dataset_image.png') 
        retrieved_features = [user_smiles, 
                                round(prediction[0][0], 3)*100]
        try:
            return render_template("pls-prediction.html", tasks = retrieved_features)
        except:
            return 'There was an error adding your task'
    if request.method == "GET":
        return render_template("pls-prediction.html")


@app.route('/publications', methods=['GET', 'POST'])
def publications():
    if request.method == 'POST':
        return render_template('publications.html')
    return render_template('publications.html')


@app.route('/contact', methods=['GET', 'POST'])
def contact():
    if request.method == 'POST':
        return render_template('contact.html')
    return render_template('contact.html')


@app.route('/literature_database', methods=['GET', 'POST'])
def literature_database():
    if request.method == 'POST':
        return render_template('literature_database.html')
    return render_template('literature_database.html')


@app.route('/legal', methods=['GET', 'POST'])
def legal():
    if request.method == 'POST':
        return render_template('legal.html')
    return render_template('legal.html')


@app.route('/what_is_osn', methods=['GET', 'POST'])
def what_is_osn():
    if request.method == 'POST':
        return render_template('what_is_osn.html')
    return render_template('what_is_osn.html')

@app.route('/plot', methods=['GET', 'POST'])
def plot():
    pass

@app.route('/enantioseparation', methods=['GET', 'POST'])
def enantioseparation():
    if request.method == 'POST':
        user_r_rejection = request.form['r_rejection']
        user_s_rejection = request.form['s_rejection']
        user_r_racemate = request.form['racemate']

        output = ""

        try:
            R1 = float(user_r_rejection)/100
            R2 = float(user_s_rejection)/100
            ratio = float(user_r_racemate)/100

            if R1 > 1 or R2 > 1:
                output = "Invalid input. Please enter rejection values equal to or lower than 100%"

            if ratio >= 1 or ratio <= 0:
                output = "Invalid input. Please enter a ratio value between 0 and 100% (excluding 0 and 100%)"
            
        except ValueError:
            output = "Invalid input. Please enter numerical values"

        if output == "":
            theta, theta2, theta3, ee_Rr, ee_Rp, ee_Sr, ee_Sp, eta_Rr, eta_Rp, eta_Sr, eta_Sp, ee2R, ee2S, eta2R, eta2S, ee3Rr, ee3Rp, ee3Sr, ee3Sp, eta3Rr, eta3Rp, eta3Sr, eta3Sp = calculate_values(user_r_rejection, user_s_rejection, user_r_racemate)
            plotting(theta, theta2, theta3, ee_Rr, ee_Rp, ee_Sr, ee_Sp, eta_Rr, eta_Rp, eta_Sr, eta_Sp, ee2R, ee2S, eta2R, eta2S, ee3Rr, ee3Rp, ee3Sr, ee3Sp, eta3Rr, eta3Rp, eta3Sr, eta3Sp)
        retrieved_features = [user_r_rejection,
                                user_s_rejection,
                                    user_r_racemate, output]
        return render_template('enantioseparation.html', tasks=retrieved_features)
    if request.method == 'GET':
        return render_template('enantioseparation.html')


if __name__ == "__main__":
    app.run(debug=False)