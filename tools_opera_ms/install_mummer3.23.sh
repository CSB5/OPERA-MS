
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
#echo $DIR
#To install mummer
echo " *** Check mummer installation";
command -v $DIR/MUMmer3.23/nucmer >/dev/null 2>&1 || {
    echo >&2 " *** Installing mummer ...";
    wget --no-check-certificate -O $DIR/download https://sourceforge.net/projects/mummer/files/mummer/3.23/MUMmer3.23.tar.gz/download;
    #mv download utils/
    cd $DIR;
    tar -xzvf download;
    cd $DIR/MUMmer3.23;
    make;
}


