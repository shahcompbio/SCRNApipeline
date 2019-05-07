FROM nceglia/base-r-scrna:latest

MAINTAINER Nicholas Ceglia <nickceglia@gmail.com>

RUN apt-get update \
  && apt-get install -y python3-pip python3-dev libffi-dev build-essential  \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
&& pip3 install --upgrade pip

RUN python3 -m pip install tensorflow
RUN python3 -m pip install PyYAML
RUN python3 -m pip install matplotlib
RUN python3 -m pip install pandas
RUN python3 -m pip install rpy2
RUN python3 -m pip install networkx
RUN python3 -m pip install dill
RUN python3 -m pip install azure
RUN mkdir /results
RUN mkdir /runs
RUN git clone https://github.com/shahcompbio/scvis.git /codebase/scvis/
RUN cd /codebase/scvis && python3 setup.py install
RUN git clone https://github.com/shahcompbio/SCRNApipeline.git /codebase/SCRNApipeline/
RUN git clone https://github.com/shahcompbio/pypeliner.git /codebase/pypeliner/
RUN cd /codebase/pypeliner && python3 setup.py install

CMD ["/bin/bash"]
