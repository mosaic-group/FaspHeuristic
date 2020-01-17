function outputData = loadHDF5(hdf5FileName)
%------------------------------------------------------------------------
% LOADHDF5 reads data from provided file name. For each group it creates 
% struct entry, and each of datasets existing in a group is put into that
% struct.
%
% Note: all datasets are treatet as a double - this is simplification but
%       a goal of this function is just to exchange numeric data easily
% 
% Author: Krzysztof Gonciarz, gonciarz@mpi-cbg.de
%------------------------------------------------------------------------
fileID = H5F.open(hdf5FileName, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
fileInfo = hdf5info(hdf5FileName);


% Start with highest level hierarchy (using only first since there is no 
% current need for more).
outputData = processLevel(fileID, fileInfo.GroupHierarchy(1), struct());

    function outString = getSuffixOfDelimiter(input, delimiter)
        idx = strfind(input, delimiter);
        if (isempty(idx))
            idx = 0;
        end
        idx = idx(end);
        outString = input((idx + 1):end);
    end

    function dataStruct = processLevel(fileId, currentHdf5Level, dataStruct)
        % First traverse all datasets on current level (in current group)
        if ~isempty(currentHdf5Level.Datasets)
            for ld = 1:length(currentHdf5Level.Datasets)
                % dataset name is the last part of path delimited with '/'
                fullDatasetPath = currentHdf5Level.Datasets(ld).Name;
                datasetName = getSuffixOfDelimiter(fullDatasetPath, '/');             
                datasetId = H5D.open(fileId, fullDatasetPath);
                try
                    
                    isString = H5T.detect_class(H5D.get_type(datasetId), 'H5T_STRING');
                    if isString == 1
                        strType = H5T.copy('H5T_C_S1');
                        H5T.set_size(strType, 256);
                        H5T.set_cset(strType, 'H5T_CSET_ASCII');
                        H5T.set_strpad(strType, 'H5T_STR_NULLTERM');
                        dataStruct.(datasetName) = cellstr(H5D.read(datasetId, strType, 'H5S_ALL','H5S_ALL','H5P_DEFAULT')');
                        H5T.close(strType);
                    else
                        % shortcut - always treat non string data as a double
                        dataStruct.(datasetName) = double(H5D.read(datasetId));
                    end
%                     dataStruct.(datasetName) = double(H5D.read(datasetId));
                    H5D.close(datasetId)
                catch
                    H5D.close(datasetId)
                end
                currentHdf5Level.Datasets(ld);
            end
        end
        
        % Go deeper if any groups available
        for level = 1:length(currentHdf5Level.Groups)
            if ~isempty(currentHdf5Level.Groups(level))                
                % path of groups starts with '/' and then each next group
                % is delimited with '.'
                fullGroupPath = currentHdf5Level.Groups(level).Name;
                groupName = getSuffixOfDelimiter(fullGroupPath, '/');
                groupName = getSuffixOfDelimiter(groupName, '.');
                
                % create new struct for a current group and dive deeper
                dataStruct.(groupName) = struct();
                dataStruct.(groupName) = processLevel(fileId, currentHdf5Level.Groups(level), dataStruct.(groupName));
            end
        end        
    end

end


