//
// Created by gonciarz on 2019-04-29.
//

#ifndef DATAHDF5_H
#define DATAHDF5_H


#include <string>
#include <map>
#include <vector>
#include <tools/prettyprint.h>
#include <hdf5.h>
#include <iostream>
#include <memory>


// Types converter from C++ to HDF5 type - won't compile if type is not defined explicitely here
template<typename DATA_TYPE> struct Hdf5Type {static hid_t type() {return  DATA_TYPE::CANNOT_DETECT_TYPE_AND_WILL_NOT_COMPILE;}};
template<> struct Hdf5Type<int8_t> {static hid_t type() {return H5T_NATIVE_INT8;}};
template<> struct Hdf5Type<uint8_t> {static hid_t type() {return H5T_NATIVE_UINT8;}};
template<> struct Hdf5Type<int16_t> {static hid_t type() {return H5T_NATIVE_INT16;}};
template<> struct Hdf5Type<uint16_t> {static hid_t type() {return H5T_NATIVE_UINT16;}};
template<> struct Hdf5Type<int> {static hid_t type() {return H5T_NATIVE_INT;}};
template<> struct Hdf5Type<unsigned int> {static hid_t type() {return H5T_NATIVE_UINT;}};
template<> struct Hdf5Type<int64_t> {static hid_t type() {return H5T_NATIVE_INT64;}};
template<> struct Hdf5Type<uint64_t> {static hid_t type() {return H5T_NATIVE_UINT64;}};
template<> struct Hdf5Type<float> {static hid_t type() {return H5T_NATIVE_FLOAT;}};
template<> struct Hdf5Type<double> {static hid_t type() {return H5T_NATIVE_DOUBLE;}};


/**
 * Easy (but good enough ;-) ) class to export data to HDF5 file format. No sophisticated checks are done - this is
 * for debug and gathering some benchmark data only - do not use it for medical and/or military purpopses.
 *
 * Example of usage:
 *
 *     DataHdf5<double> d("/tmp/out.h5");
 *     d.put("myArray", 3);
 *     d.put("sthElse", 2);
 *     d.put("sthElse", 4);
 *     d.save();
 */
template <typename T>
class DataHdf5 {
    static constexpr char DUMMY[] = "";

    const hid_t Hdf5DataType = Hdf5Type<T>::type();
    using ContainerType = std::vector<T>;

    hid_t fileId = -1;
    hid_t groupId = -1;

    std::map<std::string, ContainerType> iData;
    std::map<std::string, std::vector<std::string>> iStrData;

    void hdf5WriteData(hid_t obj_id, hid_t type_id, const char *aDataSetName, hsize_t *dims, const void *data) {

        hid_t plist_id  = H5Pcreate(H5P_DATASET_CREATE);
        hid_t space_id = H5Screate_simple(/*rank*/ 1, dims, NULL);
        hid_t dset_id = H5Dcreate2(obj_id, aDataSetName, type_id, space_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
        H5Dwrite(dset_id,type_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
        H5Dclose(dset_id);
        H5Sclose(space_id);
        H5Pclose(plist_id);
    }

    void writeString(hid_t obj_id, const char *aDataSetName,const std::vector<std::string> &s) {
        const int MaxStrLen = 256;

        // Create dataype
        hid_t dtype = H5Tcopy(H5T_C_S1);
        size_t size = MaxStrLen * sizeof(char);
        H5Tset_size(dtype, size);

        // create dataset
        hsize_t dims[] = {s.size()};
        hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
        hid_t dataset_id = H5Dcreate(obj_id, aDataSetName, dtype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // save data
        std::unique_ptr<char[]> strData(new char[MaxStrLen * s.size()]);
        memset(strData.get(), 0, MaxStrLen * s.size());
        std::size_t i = 0;
        for (auto &str : s) {
            strcpy(strData.get() + i, str.c_str());
            i += MaxStrLen;
        }
        H5Dwrite (dataset_id, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, strData.get());

        // close things
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
        H5Tclose(dtype);
    }

    template<typename CONTAINER_TYPE>
    void write(hid_t aObjectId, const std::string &aName, const CONTAINER_TYPE &aContainer) {
        hsize_t dims[] = {aContainer.size()};
        hdf5WriteData(aObjectId, Hdf5Type<typename CONTAINER_TYPE::value_type>::type(), aName.c_str(), dims, aContainer.data());
    }

public:
    explicit  DataHdf5(const std::string &aFileName, bool aDummyRun = false) {
        if (aDummyRun) return;
        fileId = H5Fcreate(aFileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (fileId == -1) {
            std::cerr << "Could not create file [" << aFileName << "]" << std::endl;
            return;
        }
        groupId = H5Gcreate2(fileId, "Analysis_data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    ~DataHdf5() {
        if (groupId != -1) H5Gclose(groupId);
        if (fileId != -1) H5Fclose(fileId);
    }

    void put(const std::string &aTableName, const T &aValue) {
        iData.try_emplace(aTableName, ContainerType{}).first->second.emplace_back(aValue);
    }

    void put(const std::string &aTableName, const std::string &aValue) {
        iStrData.try_emplace(aTableName, std::vector<std::string>{}).first->second.emplace_back(aValue);
    }

    void save() {
        if (groupId != -1) {
            for (auto& [name, data] : iData) write(groupId, name, data);
            for (auto& [name, data] : iStrData) writeString(groupId, name.c_str(), data);
         }
    }

    friend std::ostream &operator<<(std::ostream &os, const DataHdf5 &obj) {
        return os << "DataHdf5{" << "Numeric: " << obj.iData << "String: " << obj.iStrData << "}";
    }
};


#endif
