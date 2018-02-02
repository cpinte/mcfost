pushd test_data/ref3.0
rm -rf data_*
~/mcfost/src/mcfost ref3.0.para -mol
rm -f _dust_prop_th.tmp
~/mcfost/src/mcfost ref3.0.para -img 1.0
~/mcfost/src/mcfost ref3.0.para -img 10
~/mcfost/src/mcfost ref3.0.para -img 100
~/mcfost/src/mcfost ref3.0.para -img 1000
popd
