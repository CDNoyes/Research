function x = Saturate(y,low,high)

x = y;
x(x>high)= high;
x(x<low) = low;