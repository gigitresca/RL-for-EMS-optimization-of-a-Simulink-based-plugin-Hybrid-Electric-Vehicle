function generateDrivingCycleObjects(xlsx_file)
% Function to generate one DrivingCycle object to each driving cycle stored
% into the input Excel file
    in = input("Do you want to add new cycles? Y/N [Default N]: ","s");
    if isempty(in)
        in = 'N';
    end
    sheets = sheetnames(xlsx_file);
    if strcmp(in,"N")
        for i = 1:length(sheets)
            tab_data = readtable(xlsx_file,"Sheet",sheets(i),"VariableNamingRule","preserve");
            eval(sheets(i)+" = DrivingCycle(tab_data,sheets(i));");
        end
    else
        load DrivingCycleObjects.mat
        [new_cycles,~]=listdlg('PromptString','SELECT THE DRIVING CYCLES TO BE ADDED','ListString',sheets,"ListSize",[400 600],"Name",'New driving cycles selection');
        for i = new_cycles
            tab_data = readtable(xlsx_file,"Sheet",sheets(i),"VariableNamingRule","preserve");
            eval(sheets(i)+" = DrivingCycle(tab_data,sheets(i));");
        end
    end
    sheets = convertStringsToChars(sheets);
    save("data\DrivingCycleObjects",sheets{:})
end