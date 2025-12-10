FROM continuumio/miniconda3

WORKDIR /app
COPY . /app

# Create and activate environment
RUN conda update -n base -c defaults conda -y && \
    conda create -n molecul python=3.10 -y && \
    /bin/bash -lc "source $(conda info --base)/etc/profile.d/conda.sh && conda activate molecul && conda install -c conda-forge --file requirements-conda.txt -y && pip install -r requirements.txt"

EXPOSE 8501
CMD ["/bin/bash", "-lc", "source $(conda info --base)/etc/profile.d/conda.sh && conda activate molecul && streamlit run app.py --server.port=8501 --server.address=0.0.0.0"]
