#!/bin/bash
TEST_DATA_FILE="test_data_c5e7889ace157f1d440c2a5f56c8a27f7622b069_26_10_2026.tgz"
TEST_DATA_ID="1mFnJidjPK6VKrGUOrpQ21j-BP3Mdneht"
TEST_DATA_URL="https://drive.usercontent.google.com/download?export=download&id=${TEST_DATA_ID}&confirm=t"

if command -v curl &> /dev/null; then
    curl -L -o "${TEST_DATA_FILE}" "${TEST_DATA_URL}"
elif command -v wget &> /dev/null; then
    wget -nv -O "${TEST_DATA_FILE}" "${TEST_DATA_URL}"
fi

rm -rf test_data
tar xzf "${TEST_DATA_FILE}"


#Previous versions:
#test_data_7ffa5623ce83a717f3c7302f810e7ae3be3d916d_26_10_2026.tgz until 26/10/2026
#test_data_e33ad4591acbca8ab90c51a8dcf66b52e34d0d08_20_10_2026.tgz until 26/10/2026
#test_data_5580d37651d9bf2988994eb6f9cc539abc8c2561_16_10_2026.tgz until 20/10/2026
#test_data_5498a2a52a99fe085078f2c79c833d762464b0b9_21_08_2026.tgz until 16/10/2026
#test_data_0c484a1a905c13bea4c3bc0494ba288daddbe779.tgz until 21/08/2026
#test_data_b17a7fd4bce16a7c2befe713f6dd91773c057740.tar.gz
#"test_data_ca20c69767a8d3bc8da2831d258424b4ae851ab3.tgz"
