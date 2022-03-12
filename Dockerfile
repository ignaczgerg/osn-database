FROM heroku/miniconda

# Grab requirements.txt.
ADD ./requirements.txt /tmp/requirements.txt

# Install dependencies
RUN pip install -qr /tmp/requirements.txt

# Add our code
ADD ./ /opt/osn-database/
WORKDIR /opt/osn-database

RUN conda install -c rdkit rdkit pandas numpy

CMD gunicorn --bind 0.0.0.0:$PORT wsgi