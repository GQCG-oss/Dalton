def configure(options, input_files, extra_args):
    """
    This function is used by runtest to configure runtest
    at runtime for Dalton specific launch command and file naming.
    """

    from os import path

    launcher = 'dalton'
    launcher_full_path = path.normpath(path.join(options.binary_dir, launcher))

    (inp, mol) = input_files

    inp_no_suffix = path.splitext(inp)[0]
    mol_no_suffix = path.splitext(mol)[0]

    command = []

    command.append(launcher_full_path)
    command.append('-noarch -nobackup')

    if extra_args is not None:
        command.append(extra_args)

    if mol is None:
        output_prefix = inp_no_suffix
        command.append(inp_no_suffix)
    else:
        command.append('{0} {1}'.format(inp_no_suffix, mol_no_suffix))
        if inp_no_suffix == mol_no_suffix:
            output_prefix = inp_no_suffix
        else:
            output_prefix = '{0}_{1}'.format(inp_no_suffix, mol_no_suffix)

    full_command = ' '.join(command)

    relative_reference_path = 'result'

    return launcher, full_command, output_prefix, relative_reference_path
