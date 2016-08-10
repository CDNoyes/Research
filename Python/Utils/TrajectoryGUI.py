import matplotlib
matplotlib.use('TkAgg')

from numpy import arange, sin, pi
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from numpy import loadtxt

import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

root = Tk.Tk()
root.wm_title('Trajectory data visualization GUI')

# Read the available files
from os import listdir
if len(sys.argv) > 1:
    results_dir = sys.argv[1]
    if not results_dir[-1] == '/':
        results_dir = results_dir + '/'
else:
    results_dir = 'results/'
    
results_file_extension = '.txt'

results_dir_files = listdir(results_dir)
result_files = [file for file in results_dir_files if (file.find(results_file_extension)>-1 and not 'event' in file and not 'inputs' in file and not 'summary' in file)]
result_files.insert(0,result_files.pop(result_files.index('trajectory.txt'))) # move trajectory.dat to first position


# Read the headers in result_files
varlist = []
data = []
for file in result_files:
    with open(results_dir + file) as myfile:
        head = [next(myfile) for x in xrange(3)]
    
        varnames = head[1].split()
        varunits = head[2].split()
        varnameunits = [a + ', ' + b for a,b in zip(varnames,varunits)]
        varlist.append(varnameunits)
    data.append(loadtxt(results_dir+file, skiprows=3))
    

# Display dropboxes
def new_y_file(file_selected):
    file_idx = result_files.index(file_selected)
    y_var_select['menu'].delete(0, 'end') # Delete all options
    for var in varlist[file_idx]:
        y_var_select['menu'].add_command(label=var, command=Tk._setit(y_var_var,var)) # Add new options
    y_var_var.set(varlist[file_idx][0]) # Select the first one
    
    # Same with the y listbox
    y_listbox.delete(0,y_listbox.size())
    for i in range(len(varlist[file_idx])):
        y_listbox.insert(Tk.END,varlist[file_idx][i])
    y_var_last_selected = 0

def new_x_file(file_selected):
    file_idx = result_files.index(file_selected)
    x_var_select['menu'].delete(0, 'end') # Delete all options
    for var in varlist[file_idx]:
        x_var_select['menu'].add_command(label=var, command=Tk._setit(x_var_var,var)) # Add new options
    x_var_var.set(varlist[file_idx][0]) # Select the first one
    
    # Same with the x listbox
    x_listbox.delete(0,x_listbox.size())
    for i in range(len(varlist[file_idx])):
        x_listbox.insert(Tk.END,varlist[file_idx][i])
    x_var_last_selected = 0
    
def new_y_var(*args):
    var_selected = y_var_var.get()
    plot_vars(x_file_var.get(),x_var_var.get(),y_file_var.get(),y_var_var.get())
    
def new_x_var(*args):
    var_selected = x_var_var.get()
    plot_vars(x_file_var.get(),x_var_var.get(),y_file_var.get(),y_var_var.get())
    
    
def plot_vars(x_file_str,x_var_str,y_file_str,y_var_str):
    x_file_idx = result_files.index(x_file_str)
    y_file_idx = result_files.index(y_file_str)
    x_var_idx = varlist[x_file_idx].index(x_var_str)
    y_var_idx = varlist[y_file_idx].index(y_var_str)

    if y_rad2deg_var.get():
        y_multiplier = 180.0/pi
    else:
        y_multiplier = 1.0
        
    if x_rad2deg_var.get():
        x_multiplier = 180.0/pi
    else:
        x_multiplier = 1.0
        
    if y_var_str.find('[')>-1:
        # Y variable is a vector
        y_var_str_root = y_var_str.split('[')[0]
        root_found = [str.find(y_var_str_root) for str in varlist[y_file_idx]]
        idx_list=[]
        for i in range(len(root_found)):
            if not root_found[i]:
                idx_list.append(i)
                
        f.clf()
        a = f.add_subplot(111)

        x_data = data[x_file_idx] #loadtxt(results_dir+x_file_str, skiprows=3)
        y_data = data[y_file_idx] #loadtxt(results_dir+y_file_str, skiprows=3)

        colorlist =  ['r','g','b','k','y','m']*2
        legendlist = ['0','1','2','3','4','5','6','7','8','9','10','11']
        for i in range(len(idx_list)):
            a.plot(x_data[:,x_var_idx]*x_multiplier,y_data[:,idx_list[i]]*y_multiplier,colorlist[i],label=legendlist[i])
        a.grid()
        lg=a.legend()
        lg.draggable(state=True)
        canvas.show()
        toolbar.update()
        
    else:
        # Y variable is a scalar
        x_data = data[x_file_idx] #loadtxt(results_dir+x_file_str, skiprows=3)
        y_data = data[y_file_idx] #loadtxt(results_dir+y_file_str, skiprows=3)

        f.clf()
        a = f.add_subplot(111)
        a.plot(x_data[:,x_var_idx]*x_multiplier,y_data[:,y_var_idx]*y_multiplier)
        a.grid()
        canvas.show()
        toolbar.update()
    
    
y_file_var = Tk.StringVar(root)
y_file_var.set(result_files[0])
y_file_select = apply(Tk.OptionMenu,(root, y_file_var) + tuple(result_files), {'command':new_y_file})
y_file_select.pack()

y_var_var = Tk.StringVar(root)
y_var_var.set(varlist[0][0])
y_var_select = apply(Tk.OptionMenu, (root, y_var_var) + tuple(varlist[0]))
y_var_select.pack()
y_var_var.trace('w',new_y_var)
y_var_select['menu'].config(font=('calibri',(9)))


# Add y_listbox
y_frame = Tk.Frame(root, bd=2, relief=Tk.SUNKEN)
y_scrollbar = Tk.Scrollbar(y_frame)
y_scrollbar.pack(side=Tk.RIGHT, fill=Tk.Y)
y_listbox = Tk.Listbox(y_frame,height=5, width=30, selectmode='single')
for i in range(len(varlist[0])):
    y_listbox.insert(i+1,varlist[0][i])
y_listbox.pack()
y_listbox.config(yscrollcommand=y_scrollbar.set)
y_scrollbar.config(command=y_listbox.yview)
y_frame.pack()
y_var_last_selected = 0
def y_execfun(event):
    widget = event.widget
    selection=widget.curselection()
    if selection:
        # is the selection is not empty
        value = widget.get(selection[0])
        plot_vars(x_file_var.get(),x_listbox.get(x_var_last_selected),y_file_var.get(),value)
        y_var_last_selected = selection[0]
        
        # Change the dropbox value
        y_var_var.set(value)
y_listbox.activate(1)
y_listbox.bind("<Double-Button-1>",y_execfun)

# Add Y unit select checkbox
def toggle_y_rad2deg(*args):
    new_y_var()
y_rad2deg_var = Tk.IntVar()
y_rad2deg_cb = Tk.Checkbutton(root, text="RAD to DEG", variable=y_rad2deg_var,command=toggle_y_rad2deg)
y_rad2deg_cb.pack()
y_rad2deg_var.set(0)

#Display plot
f = Figure(figsize=(5,4), dpi=100)
a = f.add_subplot(111)

canvas = FigureCanvasTkAgg(f, master=root)
canvas.show()
canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

toolbar = NavigationToolbar2TkAgg( canvas, root )
toolbar.update()
canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

# X axis select
x_file_var = Tk.StringVar(root)
x_file_var.set(result_files[0])
x_file_select = apply(Tk.OptionMenu,(root, x_file_var) + tuple(result_files), {'command':new_x_file})
x_file_select.pack()

x_var_var = Tk.StringVar(root)
x_var_var.set(varlist[0][0])
x_var_select = apply(Tk.OptionMenu, (root, x_var_var) + tuple(varlist[0]))
x_var_select.pack()
x_var_var.trace('w',new_x_var)
x_var_select['menu'].config(font=('calibri',(9)))


# Add x_listbox
x_frame = Tk.Frame(root, bd=2, relief=Tk.SUNKEN)
x_scrollbar = Tk.Scrollbar(x_frame)
x_scrollbar.pack(side=Tk.RIGHT, fill=Tk.Y)
x_listbox = Tk.Listbox(x_frame,height=3, selectmode='single')
for i in range(len(varlist[0])):
    x_listbox.insert(i+1,varlist[0][i])
x_listbox.pack()
x_listbox.config(yscrollcommand=x_scrollbar.set)
x_scrollbar.config(command=x_listbox.yview)
x_frame.pack()
x_var_last_selected = 0
def x_execfun(event):
    widget = event.widget
    selection=widget.curselection()
    if selection:
        # is the selection is not empty
        value = widget.get(selection[0])
        plot_vars(x_file_var.get(),value,y_file_var.get(),x_listbox.get(y_var_last_selected))
        x_var_last_selected = selection[0]
        
        # Change the dropbox value
        x_var_var.set(value)

x_listbox.activate(1)
x_listbox.bind("<Double-Button-1>",x_execfun)

# Add X unit select checkbox
def toggle_x_rad2deg(*args):
    new_x_var()
x_rad2deg_var = Tk.IntVar()
x_rad2deg_cb = Tk.Checkbutton(root, text="RAD to DEG", variable=x_rad2deg_var,command=toggle_x_rad2deg)
x_rad2deg_cb.pack()
x_rad2deg_var.set(0)



# Execute once to load variable list correctly
new_y_file(result_files[0]) 
new_x_file(result_files[0])

# Run!
Tk.mainloop()

