#!/bin/bash
TEST_DATA_FILE="test_data_5498a2a52a99fe085078f2c79c833d762464b0b9_21_08_2026.tgz"
TEST_DATA_ID="1p8rml0kkkarFNDN7no0SDrr-kWvJ_QID"
TEST_DATA_URL="https://drive.usercontent.google.com/download?export=download&id=${TEST_DATA_ID}&confirm=t"

if command -v curl &> /dev/null; then
    curl -L -o "${TEST_DATA_FILE}" "${TEST_DATA_URL}"
elif command -v wget &> /dev/null; then
    wget -nv -O "${TEST_DATA_FILE}" "${TEST_DATA_URL}"
fi

rm -rf test_data
tar xzf "${TEST_DATA_FILE}"


#Previous versions:
#test_data_0c484a1a905c13bea4c3bc0494ba288daddbe779.tgz until 21/08/2026
#test_data_b17a7fd4bce16a7c2befe713f6dd91773c057740.tar.gz
#"test_data_ca20c69767a8d3bc8da2831d258424b4ae851ab3.tgz"
