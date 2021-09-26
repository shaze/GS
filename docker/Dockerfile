FROM quay.io/cybozu/ubuntu:20.04

ENV PERL_MM_USE_DEFAULT=1

RUN apt-get -y update && \
    echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections  &&\
    apt-get -y install perl make gcc

RUN apt-get -y update && apt-get -y install expat emacs-nox libexpat1-dev libxml-sax-expat-perl

ENV PERL5LIB=/usr/share/modules/init

RUN cpan -i XML::Parser::PerlSAX XML::DOM XML::Twig &&\
    cpan -i Encode::Locale File::Listing HTML::Entities HTML::HeadParser HTTP::Cookies HTTP::Date HTTP::Negotiate HTTP::Request HTTP::Request::Common HTTP::Response HTTP::Status LWP::MediaTypes Net::HTTP Test::Fatal Test::RequiresInternet URI::Escape WWW::RobotRules Test::Needs Try::Tiny LWP::MediaTypes URI URI::URL HTML::Tagset HTTP::Headers HTTP::Date IO::HTML && \
    cpan -i Term::ReadLine
    


RUN apt-get -y update && \
    apt-get install -y apt-utils  &&\
    apt-get install -y locales wget unzip python3.9-full python2.7 perl git && \
    ln -s /usr/bin/python2.7 /usr/bin/python2 && \
    rm -rf /var/lib/apt/lists/* \
    && localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8



ENV LANG en_US.utf8


RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip\
    && unzip fastqc_v0.11.9.zip && \
    mv FastQC /usr/local && ln -s /usr/local/FastQC/fastqc /usr/local/bin/fastqc && \
    chmod a+x /usr/local/bin/fastqc

RUN wget https://bootstrap.pypa.io/get-pip.py && python3.9 get-pip.py &&\
    /bin/rm get-pip.py && \
    pip3.9 install multiqc==1.11 cutadapt==3.4 pandas &&\
    wget https://bootstrap.pypa.io/pip/2.7/get-pip.py &&\
    python2 get-pip.py  
    

RUN git clone https://github.com/katholt/srst2 &&\
    cd srst2 && git checkout v0.2.0 && cd .. && pip2.7 install srst2/ &&\
    pip2.7 install scipy

RUN wget https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux &&\
    mv prodigal.linux /usr/local/bin/prodigal && chmod a+x /usr/local/bin/prodigal

RUN  wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.12.0+-x64-linux.tar.gz &&\
     tar -xzf ncbi-blast-2.12.0+-x64-linux.tar.gz && \
     mv  ncbi-blast-2.12.0+/bin/* /usr/local/bin

RUN  wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary &&\
     mv bedtools.static.binary /usr/local/bin/bedtools &&\
     chmod a+x /usr/local/bin/bedtools


RUN  wget https://github.com/freebayes/freebayes/releases/download/v1.3.4/freebayes-1.3.4-linux-static-AMD64.gz &&\
     mv freebayes-1.3.4-linux-static-AMD64.gz /usr/local/bin/freebayes && \
     chmod a+x /usr/local/bin/freebayes 




RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/bowtie2-2.2.9-linux-x86_64.zip && \
    unzip bowtie2-2.2.9-linux-x86_64.zip && \
    mv bowtie2-2.2.9/bowtie2* /usr/local/bin/

RUN apt-get -y update && apt-get install -y samtools=1.10-3 velvet=1.2.10+dfsg1-7

RUN apt-get -y install environment-modules=4.4.1-1 &&\
    ln -s /usr/bin/python2.7 /usr/bin/python

RUN apt-get -y autoremove  & /bin/rm -rf /root/*


     
