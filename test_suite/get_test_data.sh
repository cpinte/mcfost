TEST_DATA_FILE="test_data_ca20c69767a8d3bc8da2831d258424b4ae851ab3.tgz"
TEST_DATA_URL="https://drive.usercontent.google.com/download?export=download&id=1_slLnqLifQJzHpTrpbuCYVyUoGFpDqoZ&confirm=t"

wget -O "${TEST_DATA_FILE}" "${TEST_DATA_URL}"
rm -rf test_data
tar xzf "${TEST_DATA_FILE}"

#Previous versions:
#test_data_0c484a1a905c13bea4c3bc0494ba288daddbe779.tgz
#test_data_b17a7fd4bce16a7c2befe713f6dd91773c057740.tar.gz
