%Start workers on other hosts
clientHost = 'mipttea';
node = {'mipttea', 'teawin7'};
for i = 1:length(node)
    for j = 1:8
        str = ['!startworker -name w_' num2str(j) '_' node{i} ' -jobmanagerhost ' clientHost ' -jobmanager jm -remotehost ' node{i} ' -v'];
        eval(str)
    end
end