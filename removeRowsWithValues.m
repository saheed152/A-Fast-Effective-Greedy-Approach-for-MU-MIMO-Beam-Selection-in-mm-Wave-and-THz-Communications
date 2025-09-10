function [modifiedMatrix] = removeRowsWithValues(matrix, valuesToRemove)
    % This function removes the rows from 'matrix' which contain any element of 'valuesToRemove'

    % Argument 'matrix' is a two-dimensional matrix.
    % Argument 'valuesToRemove' is a vector containing the values that, if present in a row, 
    % will lead to the deletion of the entire row from 'matrix'.
    
    % Initialize a logical index vector to mark rows for deletion
    rowsToDelete = false(size(matrix, 1), 1);
    
    % Loop through each value in valuesToRemove
    for value = valuesToRemove
        % Find rows that contain the current value
        rowsWithCurrentValue = any(matrix == value, 2);
        
        % Combine with previous results
        rowsToDelete = rowsToDelete | rowsWithCurrentValue;
    end
    
    % Delete the marked rows
    matrix(rowsToDelete, :) = [];
    
    % Return the modified matrix
    modifiedMatrix = matrix;
end
