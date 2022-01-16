clear;
clear classes;

if count(py.sys.path,'') == 0
    insert(py.sys.path,int32(0),'');
end

mod = py.importlib.import_module('mymodule')
py.reload(mod)

py.mymodule.add_numbers(2, 2)
py.mymodule.add_list({1, 2, 3})
py.mymodule.add_list(1:5)