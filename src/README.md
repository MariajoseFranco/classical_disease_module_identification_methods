# Disease Module Identification - Classical Methods
**By:** Mariajose Fraco Orozco

**Tutor:** Lucia Prieto Santamaria

To run it locally first create a virtual environment and install all the dependencies required.

## Create a Virtual Environment
### In Mac:
```bash
python3 -m venv <myenvname>
```
Then, to activate it:
```bash
source <myenvname>/bin/activate
```
and to deactivate it
```bash
deactivate
```
### In Windows:
```bash
python -m venv <myenvname>
```
Then, to activate it:
```bash
<myenvname>\Scripts\activate.bat
```
and to deactivate it
```bash
deactivate
```
## Install the dependencies

Once the virtual environment is activated, install the dependencies with the following command
```bash
pip install -r requirements.txt
```

## Run main.py
When everything is set up, you can run the different algorithms locally or in docker.
### For running it locally...
```bash
python main.py
```
### For running it in docker desktop...
To create a docker image:
```bash
docker build -t <project_name> .
```

To run the docker image created:
```bash
docker run -v "$(pwd)/outputs":/app/outputs <project_name>
```