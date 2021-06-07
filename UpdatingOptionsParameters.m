function[options_updated,n_changed_fields,changed_fields,stat]= UpdatingOptionsParameters( options_updated, options_argin, updatetype, isdebug )
    if nargin<1
          options_updated.q=[1:5]';
          options_updated.ww='ddfsdfs';
          options_updated.qqq=1111;
          options_argin.q = fieldnames(options_updated);   
          options_argin.q22 = fieldnames(options_updated);   
          updatetype =  'exact'   ;    
    end
    stat.n_changed_fields = []; 
    stat.n_updated_fields = []; 
    stat.n_added_fields = []; 
    stat.n_discarded_fields= []; 
    stat.n_deleted_fields = []; 
    stat.added_fields = []; 
    stat.updated_fields = []; 
    stat.changed_fields = [];       
    stat.discarded_fields = [];     
    n_changed_fields = 0;
    changed_fields  = [];    
    if ~isa(options_updated,'struct'); error(['options_updated must be a struct variable'] );  end 
    if ~exist('options_argin','var')||isempty(options_argin)
        AllFieldNameSet   = fieldnames( options_updated ); 
        FieldNameSet_argin   = []; 
        FieldNameSet_updated = []; 
        n_changed_fields     = 0; 
        disp('There is no updating of options.');
        return; 
    end
    if ~isa(options_argin,'struct') ; error(['options_argin must be a struct variable'] );  end 
    if ~exist('updatetype','var')||isempty(updatetype)  
        updatetype = 'exact'; 
    end
    if ~exist('isdebug','var')||isempty(isdebug)  
        isdebug = true; 
    end
    AllFieldNameSet      = fieldnames( options_updated );  
    FieldNameSet_argin   = fieldnames( options_argin  );  % these fields will be updated. 
    FieldNameSet_updated = intersect( AllFieldNameSet ,  FieldNameSet_argin ); 
    otherfieldnames      = setdiff(FieldNameSet_argin, AllFieldNameSet); 
    switch lower(updatetype)
        case 'exact' 
            if ~isempty( otherfieldnames )
                disp( otherfieldnames  )
                error( ['The above field(s) does not exsit in fieldnames of options. The updating is failing.'] );
            end 
            n_updated_fields = length(FieldNameSet_updated); 
            for i = 1: n_updated_fields
                sname = FieldNameSet_updated{i};  
                options_updated.(sname) = options_argin.(sname);   
            end 
            stat.n_updated_fields = n_updated_fields; 
            stat.n_changed_fields = n_updated_fields; 
            stat.updated_fields   = FieldNameSet_updated; 
            stat.changed_fields   = FieldNameSet_updated; 
            if isdebug 
                if n_updated_fields>0
                    disp('Update the following fields of options:');%  disp( FieldNameSet_updated );                
                    [options_temp]=extractVariablesFromFieldnames(options_updated, FieldNameSet_updated ) ;
                    disp( options_temp ); 
                else
                    disp( 'No updating of fields.' );
                end 
            end 
        case {'discard', 'ignore'}
            n_discarded_fields = length( otherfieldnames  ); 
            n_updated_fields = length(FieldNameSet_updated); 
            for i = 1: n_updated_fields
                sname = FieldNameSet_updated{i};  
                options_updated.(sname) = options_argin.(sname);   
            end 
            stat.n_changed_fields  = n_updated_fields; 
            stat.n_updated_fields  = n_updated_fields; 
            stat.n_discarded_fields = n_discarded_fields; 
            stat.updated_fields   = FieldNameSet_updated; 
            stat.changed_fields   = FieldNameSet_updated; 
            stat.discarded_fields = otherfieldnames; 
            if isdebug 
                if n_updated_fields>0
                    disp( 'Update the following fields of options:'  ); %  disp( FieldNameSet_updated );                
                    [options_temp] = extractVariablesFromFieldnames(   options_updated, FieldNameSet_updated )  ;
                    disp( options_temp ); 
                else
                    disp( 'No updating of fields.'  );
                end
                if n_discarded_fields>0 
                    disp('Discard/ignore the following fields of options_argin:' );      
                    [options_temp]=extractVariablesFromFieldnames(options_argin, otherfieldnames);
                    disp( options_temp ); 
                else
                    disp( 'No adding of fields.'  );
                end
            end 
        case 'add' 
            n_added_fields = length( otherfieldnames  ); 
            for i = 1: n_added_fields
                sname = otherfieldnames{i};  
                options_updated.(sname) = options_argin.(sname);   
            end   
            n_updated_fields     = length(FieldNameSet_updated); 
            for i = 1: n_updated_fields
                sname = FieldNameSet_updated{i};  
                options_updated.(sname) = options_argin.(sname);   
            end 
            n_changed_fields = n_updated_fields + n_added_fields ; 
            stat.n_added_fields = n_added_fields; 
            stat.n_updated_fields = n_updated_fields; 
            stat.n_changed_fields = n_changed_fields; 
            stat.added_fields     = otherfieldnames; 
            stat.updated_fields   = FieldNameSet_updated; 
            stat.changed_fields   = FieldNameSet_argin;  
            if isdebug 
                if n_updated_fields>0
                    disp( 'Update the following fields of options:'  );               
                    [options_temp] = extractVariablesFromFieldnames(options_updated, FieldNameSet_updated ) ;
                    disp( options_temp ); 
                else
                    disp(  'No updating of fields.'  );
                end
                if n_added_fields>0 
                    disp( 'Add the following fields into options:'  );              
                    [options_temp] = extractVariablesFromFieldnames(   options_updated, otherfieldnames )  ;
                    disp( options_temp ); 
                else
                    disp(  'No adding of fields.'  );
                end
            end 
        case 'delete'  
            n_deleted_fields = length( FieldNameSet_updated ); 
            if ~isempty( FieldNameSet_updated ) && isdebug 
                disp( FieldNameSet_updated  );  
                warning(  ['delete the above field(s) that exsit in fields of options.']  );
            end  
            options_updated = rmfield(options_updated, FieldNameSet_updated );   
            n_changed_fields = n_deleted_fields; 
            stat.n_deleted_fields = n_deleted_fields; 
            stat.n_changed_fields = n_deleted_fields; 
            stat.n_updated_fields = []; 
            stat.n_added_fields   = []; 
            stat.deleted_fields   = FieldNameSet_updated; 
            stat.changed_fields   = FieldNameSet_updated;             
        otherwise
            error(  ['There is no definition of updatetype:',updatetype]  );
    end
    changed_fields = stat.changed_fields; 
end
