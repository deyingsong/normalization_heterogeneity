function merged = merge_structs(struct1, struct2)
% MERGE_STRUCTS - Combine two structures
%
% H1 Line: Merge fields from two structures
%
% Syntax:
%   merged = merge_structs(struct1, struct2)
%
% Description:
%   Combines fields from two structures, with struct2 fields
%   overwriting struct1 fields if duplicates exist.
%
% Inputs:
%   struct1 - First structure
%   struct2 - Second structure
%
% Outputs:
%   merged - Combined structure
%
% Example:
%   merged = merge_structs(params1, params2);

    merged = struct1;
    fields = fieldnames(struct2);
    for i = 1:length(fields)
        merged.(fields{i}) = struct2.(fields{i});
    end
end
