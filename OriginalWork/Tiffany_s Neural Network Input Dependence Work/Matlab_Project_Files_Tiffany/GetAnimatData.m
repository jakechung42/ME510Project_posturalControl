%Reads in a sentence of data from Animatlab. Data format is such that:
%xFF xFF - Begin message
%x01 - Message type is data
%xXX xXX - Message size 

%xXX xXX - Data ID
%xXX xXX xXX xXX - Data
%Repeat for all data

%xXX - Checksum

function newValues = GetAnimatData(r,oldValues,id_list,currentCommand)

%new values will be the old values with any new incoming data overwriting
%the old data
newValues = oldValues;
chksm = 0;

%Confirm start of data packet
for ii = 1:2 %run this code twice
    temp = fread(r,1);
    if temp ~= 255; %If we do not get xFF twice in a row kick out.
        return
    end
    chksm = chksm + temp;
end

%Confirm we are getting a packet of data type
temp = fread(r,1);
if temp ~= 1; %If we do not get 1 kick out.
    return
end
chksm = chksm + temp;

%Determine the number of data points we are getting
temp2 = fread(r,2);
message_length = typecast(uint8(temp2),'UINT16');
num_datapoints = message_length/6-1; %There are six bytes per data point (2 ID, 4 Data) and 6 bytes of non-data

ID = zeros(1,num_datapoints);
val = zeros(1,num_datapoints);

chksm = chksm + sum(temp2);

%Read in all the data
for jj = 1:num_datapoints
    
    temp3 = fread(r,2);
    ID(jj) = typecast(uint8(temp3),'UINT16');
    temp4 = fread(r,4);
    val(jj) = typecast(uint8(temp4),'single');
    
    chksm = chksm + sum(temp3);
    chksm = chksm + sum(temp4);

end

%Perform checksum check
Checksum = fread(r,1);
if Checksum ~= mod(chksm,256) % checksum
    warning('checksum error! cmd %d',currentCommand)
    return
end

%If successful checksum, write the new data points into vector
for kk = 1:num_datapoints
    j = id_list == ID(kk); %Find the index associated with the current ID value
    newValues(j) = val(kk); %Write the data into the correct index
end

end %End Program