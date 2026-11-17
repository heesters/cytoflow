from PyInstaller.utils.hooks import collect_data_files

# python-pptx starts every presentation from the template files that ship
# inside the package.
datas = collect_data_files('pptx')
