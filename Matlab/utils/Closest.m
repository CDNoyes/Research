function [val,idx] = Closest(array,value)

[~,idx] = min( abs(array-value) ); 
val = array(idx);

end