# Python script to create a quick and dirty classification web page
# Generates a few html + php files which can then be transferred to
# a web page (which runs php).
# Eric Emsellem  - ESO/CRAL 2020
# based and expanded from an idl version from Jesse van de Sande/Luca Cortese

import glob
import os
import re
from os.path import join as joinpath


# version 0.4 = Added brief doc and a few more options
# version 0.3 = Added person_name. Tested on real case
# version 0.2 = tested that it runs with now more than 1 batch.
# version 0.1 = tested that it runs. Still need to test the php / web page

def make_index_html(folder, list_classes, list_explanations, 
                    title, ntargets, nbatches, email, person_name):
    f = open(joinpath(folder, "index.html"), "w+")    
    htmltext = (f"<html><head>\n"
                f"<meta http-equiv=\"content-type\" content=\"text/html; charset=UTF-8\"><title>{title}</title></head> \n"
                f"\n"
                f"<body>\n" 
                f"<h2>Material for the classification of targets</h2>\n"
                f"\n"
                f"<p> There are {len(list_classes)} categories we wish to separate:\n"
                f"<ol>\n")

    for i, (button, expl) in enumerate(zip(list_classes, list_explanations)):
        htmltext += (f"<li><b> Class {i}: {button} = {expl} </b> </li>\n")
    htmltext += (f"</ol> \n"
                 f"</p>\n"
                 f"\n"
                 f"<br>\n"
                 f"The total sample includes {ntargets} galaxies, divided into {nbatches} batches: <br>\n"
                 f"---------<br>\n")
    for i in range(nbatches):
        htmltext += (f"<br>\n"
                     f"<b><a href=\"batch_{i+1:02d}.html\"> Batch {i+1}  </a> <br>\n"
                     f"\n"
                     f"  \n")
    htmltext += (f"---------<br>  \n"
                 f"<br>\n"
                 f"Please, complete all four batches and send the results to "
                 f"<a href=\"mailto:{email}\">{person_name}</a><br>\n"
                 f"Please send you results as ascii file (not embedded into an email) "
                 f"named as YourInitials_batchNN.txt (e.g., ee_batch1.txt)\n"
                 f"----------<br>\n"
                 f"<br>\n"
                 f"<br>\n"
                 f" \n"
                 f"</body></html>\n")
    f.write(htmltext)
    f.close()

def get_header(title, subtitle):
    return (f'<html> \n'
            f'<head><title>{title}</title></head> \n'
            f'\n'
            f'<body> \n'
            f'<h2>{subtitle}</h2> \n'
            f'\n'
            f'<form action="php/form.php" method="post" target="_blank"> \n'
            f'<TABLE cellspacing=3 cellpadding=3>\n')

def get_footer(list_classes=[]):
    footer = (f'</TABLE><br>\n'
              f'Once you have classified all galaxies in the list, '
              f'please check your work by clicking on the {len(list_classes)} buttons below <br>\n'
              f'You will be able to compare all galaxies that you have assigned '
              f'to the same type and make sure that your classification is correct. <br><br>\n'
              f'\n'
              f'\n'
              f'<br>\n')
    for button in list_classes:
        footer += (f'<input type=\"submit\" value=\"{button}\" '
                   f'name=\"submit_button\">\n<br>\n')
    footer += ('<br>\nHave you checked your classification using the '
               'above links?  <br> '
               'Are you happy with it? <br><br>\n'
               'If so, please click the button below to print your results. <br><br>\n'
               '<input type=\"submit\" value=\"Print\" name=\"submit_button\"></form> \n'
               '</body> \n'
               '</html> \n')
    return footer

def get_figtext(fig_name, target_id="dummy", numfig=0, 
                nfig_per_line=1, list_classes=[], height=350):
    if nfig_per_line == 2:
        if numfig%2 == 0: 
            text = (f'<tr><td><b>{target_id}</b></br><input type=\"hidden\" '
                    f'name=\"name[]\" value=\"+{target_id}\">\n')
        else: 
            text = (f'<td><b>{target_id}</b></br><input type=\"hidden\" '
                    f'name=\"name[]\" value=\"{target_id}\">\n')
    else:
        text = (f'<tr><td><b>{target_id}</b></br><input type=\"hidden\" '
                f'name=\"name[]\" value=\"{target_id}\">\n')

    for i, button in enumerate(list_classes):
        text += (f'<input type=\"radio\" name=\"type[{numfig}]\" alt=\"radio\" '
                 f'value=\"{i}\"> {button} &nbsp&nbsp \n')
    text += (f'<input type=\"hidden\" name=\"img[]\" '
             f'value=\"../{fig_name}\">\n' 
             f'<img height={height} src=\"{fig_name}\"> </td>\n</tr>\n')
    return text

def safely_create_folder(path, verbose=True):
    """Create a folder given by the input path
    This small function tries to create it and if it fails
    it checks whether the reason is because it is not a path
    and then warn the user
    and then warn the user
    """
    if path is None :
        if verbose : print_info("Input path is None, not doing anything")
        return
    if verbose : 
        print(f"Trying to create {path} folder", end='')
    try: 
        os.makedirs(path)
        if verbose:
            print("... Done\n", end='\n')
    except OSError:
        if not os.path.isdir(path):
            print("ERROR: Failed to create folder! Please check the path")
            return
        if os.path.isdir(path):
            if verbose:
                print("... Folder already exists, doing nothing.")

def make_form_php(folder, list_classes):
    safely_create_folder(joinpath(folder, "php/"))
    f = open(joinpath(folder, "php/form.php"), "w+")
    phptext = (f"<?php\n"
               f"\n"
               f"\n"
               f"\n"
               f"if(count($_POST[\'name\'])!=count($_POST[\'type\'])){{\n"
               f"\n"
               f"  echo \"<p> <b> Error: You have not classified all the "
               f"objects in the list! <br><br>\";\n"
               f"  echo \"You have classified \".count($_POST[\'type\']).\" "
               f"out of \".count($_POST[\'name\']).\" galaxies <br>\";\n"
               f"                    echo \"Please go back and make sure that "
               f"each galaxy has been classified. <br>\n"
               f"                        Thanks for your collaboration. </b></p>\" ;\n"
               f"                        }}\n"
               f"\n"
               f"else {{\n"
               f"\n"
               f"$name = $_POST[\'name\'];\n"
               f"$type = $_POST[\'type\'];\n"
               f"$img = $_POST[\'img\'];\n"
               f"\n")

    for i, name in enumerate(list_classes):
        namel = name.lower()
        if i == 0:
            phptext += (f"if($_POST[\'submit_button\']==\'{name}\'){{\n"
                        f"                                include \'{namel}.php\';\n"
                        f"                                }}\n"
                        f"\n")
        else:
           phptext += (f"else if($_POST[\'submit_button\']==\'{name}\'){{\n"
                       f"                                include \'{namel}.php\';\n"
                       f"                                }} \n")
    phptext += (f"else if($_POST[\'submit_button\']==\'Print\'){{\n"
                f"                                include \'print.php\';\n"
                f"                                }} \n"
                f"\n"
                f"}}\n"
                f" ?>\n")
    f.write(phptext)
    f.close()

def make_info_php(folder):
    safely_create_folder(joinpath(folder, "php/"))
    f = open(joinpath(folder, "php/info.php"), "w+")
    f.write("<?php phpinfo();?>\n")
    f.close()
        
def make_print_php(folder):
    safely_create_folder(joinpath(folder, "php/"))
    f = open(joinpath(folder, "php/print.php"), "w+")
    phptext = ("<?php \n"
               "\n"
               "\n"
               " for ($i=0;$i< count($name);$i++) {\n"
               "\n"
               "echo \" \".$name[$i].\" \".$type[$i].\"<br>\";"
               "\n"
               "}\n"
               "\n"
               " ?>\n")
    f.write(phptext)
    f.close()

def make_class_php(folder, classname, class_number, width=550):
    safely_create_folder(joinpath(folder, "php/"))
    f = open(joinpath(folder, f"php/{classname.lower()}.php"), "w+")
    phptext = (f"<html> \n"
               f"<head><title>{classname}</title></head>\n"
               f"\n"
               f"<body>\n"
               f"<h2> {classname} </h2>\n"
               f"\n"
               f"<TABLE cellspacing=3 cellpadding=3>\n"
               f"<tr>\n"
               f"\n"
               f"<?php \n"
               f"\n"
               f"\n"
               f"$count=1;\n"
               f"\n"
               f" for ($i=0;$i<count($name);$i++) {{\n"
               f"\n"
               f"        if($type[$i]=={class_number}) {{\n"
               f"\n"
               f"                 echo \"<td>\";  \n"
               f"                 echo \"<b>\".$name[$i].\"</b></br>\";\n"
               f"                 echo \"<img width={width} src=\".$img[$i].\"></td>\";\n"
               f"                 $count++;\n"
               f"                }}\n"
               f"        if ($count==3){{\n"
               f"                echo \"</tr>\";\n"
               f"                $count=1;\n"
               f"                }}\n"
               f"}}\n"
               f"\n"
               f" ?>\n"
               f"\n"
               f"</body>\n"
               f"</html>\n")
    f.write(phptext)
    f.close()

def goto_folder(newpath):
    try:
        prev_folder = os.getcwd()
        newpath = os.path.normpath(newpath)
        os.chdir(newpath)
        print(f"Going to folder {newpath}")
    except OSError:
        if not os.path.isdir(newpath):
            raise
    return prev_folder, newpath

# Main routine
def make_classification_page(folder='/misc/HomePage/Class_Maker/',
                             fig_folder="Figures/",
                             filename="batch", nbatches=1,
                             extension="png",
                             prefix="",  suffix="", nfig_per_line=1,
                             list_classes=['class1', 'class2', 'class3'],
                             list_explanations=['this class', "this other class", "this late class"],
                             title="Classification of Targets", subtitle="My subtitle",
                             email="eric.emsellem@eso.org", person_name="Eric",
                             height_per_fig=350, width_per_fig=550):
    """Generate the html and php files

    Input
    -----
    folder (str): root folder where the script will be created
    fig_folder (str): folder where the figures are stored, subfolder of "folder"
    filename (str): base name for the html files
    prefix (str): prefix used to identify the figures
    suffix (str): suffix used to identify the figures
    extension (str): extension for the figures (e.g., 'png')
    nbatches (int): number of batches (all figures will be split in nbatches)
    list_classes (list of str): names of the classes - make these explicit!
    list_explanations (list of str): explanations for the classes
    title (str): title for the page
    subtitle (str): subtitle for the page
    email (str): email address - who we should send the result to
    person_name (str): name of the person we should send this to
    height_per_fig (int): height in pixel for each figure in the batch
    width_per_fig (int): width in pixel for each figure in the summary per class

    Creates
    -------
    an index.html file which should be where to point users towards.
    batch html files containing the list of figures and a form to classify them.
    php folder with one php script per class and a few others (info, form, print).

    """
    
    # Test folder existence
    if not os.path.isdir(folder):
        print(f"ERROR: input folder {folder} does not exist")
        return 

    figures_folder = joinpath(folder, fig_folder)
    list_figures = glob.glob(f"{figures_folder}/*png")
    
    # Detect which files have the right suffix/prefix/extension
    ids = []
    for name in list_figures:
        path, fname = os.path.split(name)
        subname = re.findall(f'{prefix}(.*?){suffix}\.{extension}', fname)
        if len(subname) == 0:
            continue
        ids.append(subname[0])

    # Number of targets to classify
    ntargets = len(ids)
    print(f"Found {ntargets} targets to classify")
    
    # start each batch in a loop
    start_batch = 0
    ntargets_per_batch = ntargets // nbatches
    end_batch = ntargets_per_batch
    print("Starting writing up the batch html files", end='')
    for ib in range(nbatches):
        # Opening the html batch file
        namebatch = f"{filename}_{ib+1:02d}.html"
        f = open(joinpath(folder, namebatch), "w+")
        # Header
        f.write(get_header(title, f"{subtitle} - batch{ib+1:02d}"))
        # loop on figures
        count = 0
        for n in range(start_batch, end_batch):
            # Id of the figure
            this_id = ids[n]
            # Name of the figure
            name_figure = f"{prefix}{this_id}{suffix}.{extension}"
            # Test existence
            full_name_figure = joinpath(figures_folder, name_figure)
            if not os.path.isfile(full_name_figure):
                print(f"WARNING: file {name_figure} does not exist in folder {figures_folder}")
                continue
            # Add figure reference to html
            newtext = get_figtext(joinpath(fig_folder,name_figure), this_id,
                                 count, nfig_per_line=nfig_per_line, 
                                 list_classes=list_classes,
                                 height=height_per_fig)
            f.write(newtext)
            count += 1

        f.write(get_footer(list_classes))
        # closing the batch file
        f.close()
        # Change the range for the next batch
        end_batch += ntargets_per_batch
        if end_batch > ntargets: 
            end_batch = ntargets
        start_batch += ntargets_per_batch
        if start_batch >= ntargets:
            break
    print("...Done")

    # make the buttons php
    print("Writing up the php for each class")
    for i, button in enumerate(list_classes):
        make_class_php(folder, button, i, width=width_per_fig)

    # make the other php (form, info, print, index.html)
    print("Writing up the print php")
    make_print_php(folder)
    print("Writing up the info php")
    make_info_php(folder)
    print("Writing up the form php")
    make_form_php(folder, list_classes)
    print("Finally... Writing up the index.html file")
    make_index_html(folder, list_classes, list_explanations, 
                    title, ntargets, nbatches, email, person_name)
