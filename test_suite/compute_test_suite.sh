export mcfost=$HOME"/mcfost/src/mcfost"

for dir in test_data/*; do
    export param=`basename "$dir".para`
    pushd "$dir"
    rm -rf data_*
    $mcfost $param -mol
    rm -f _dust_prop_th.tmp
    $mcfost $param -img 1.0
    $mcfost $param -img 10
    $mcfost $param -img 100
    $mcfost $param -img 1000
    popd
done
