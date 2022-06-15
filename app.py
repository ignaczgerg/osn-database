from flask import Flask, render_template, url_for, request, redirect
from mordred._base.descriptor import Descriptor
from similarity_calc import DataLoader, Similarity
import predictor
from rdkit import Chem
from rdkit.Chem import Draw
from thetaplots import calculate_diagrams
import io
import random
from flask import Response
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

app = Flask(__name__)


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


@app.route('/enantiomer_separation', methods=['GET', 'POST'])
def plot_png():
    fig = create_figure()
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='static/png')

def create_figure():
    fig = Figure()
    axis = fig.add_subplot(1, 1, 1)
    xs = range(100)
    ys = [random.randint(1, 50) for x in xs]
    axis.plot(xs, ys)
    return fig

def enantiomer_separation():
    if request.method == 'POST':
        user_r_rejection = request.form['r_rejection']
        user_s_rejection = request.form['s_rejection']
        user_r_racemate = request.form['racemate']
        try:
            state = calculate_diagrams(user_r_rejection, user_s_rejection, user_r_racemate)
        except:
            state = False
        retrieved_features = [user_r_rejection,
                                user_s_rejection,
                                user_r_racemate]

        return render_template('enantiomer_separation.html', tasks=retrieved_features)
    return render_template('enantiomer_separation.html')


# def pls_prediction():
#     if request.method == 'POST':
#         user_smiles = request.form['content']
#         #comparing_dataset = DataLoader.csv_loader('data/osn_data_similarity.csv')
#         try:
#         #    similar = Similarity.rejection_similarity(comparing_dataset, reference=user_smiles)
#             user_descr = predictor.descripter(user_smiles)
#             prediction = predictor.predictor(user_descr)
#         except:
#             return 'Wrong SMILES, please try again. Currently, molecules with more than 5 heavy atoms work only.'

#         Draw.MolToFile(Chem.MolFromSmiles(user_smiles),'static/dataset_image.png') 
#         retrieved_features = [user_smiles, 
#                                 round(prediction[0][0], 3)*100]
#         try:
#             return render_template("pls-prediction.html", tasks = retrieved_features)
#         except:
#             return 'There was an error adding your task'
#     if request.method == "GET":
#         return render_template("pls-prediction.html")


if __name__ == "__main__":
    app.run(debug=True)