
import os
import sys
import runtest_v1 as runtest


def write_stderr(log_file, s):
    """
    Writes s to stderr and to file log_file
    unless log_file is None.
    """
    if log_file:
        with open(log_file, 'w') as f:
            f.write(s)
    sys.stderr.write(s)


class Filter(runtest.Filter):

    def __init__(self):
        runtest.Filter.__init__(self)

    def add(self, *args, **kwargs):
        try:
            runtest.Filter.add(self, *args, **kwargs)
        except runtest.FilterKeywordError as e:
            sys.stderr.write(str(e))  # FIXME currently not written to any log file
            sys.exit(-1)


class TestRun(runtest.TestRun):

    def __init__(self, _file, argv):
        runtest.TestRun.__init__(self, _file, argv)
        self.return_code = 0

    def run(self, inp_files, mol_files=[], other_files=[], f=None, args='', accepted_errors=[], noarch=True):

        launch_script = os.path.normpath(os.path.join(self.binary_dir, 'dalton'))
        if self.skip_run:
            sys.stdout.write('\nskipping actual run\n')
        else:
            if not os.path.exists(launch_script):
                sys.stderr.write('ERROR: launch script %s not found\n' % launch_script)
                sys.stderr.write('       have you set the correct --binary-dir (or -b)?\n')
                sys.stderr.write('       try also --help\n')
                sys.exit(-1)

        if noarch:
            launcher = '%s -ow -noarch -nobackup %s' % (launch_script, args)
        else:
            launcher = '%s -ow -nobackup %s' % (launch_script, args)

        commands = []
        messages = []
        outputs_no_suffix = []
        for inp in inp_files:
            inp_no_suffix = os.path.splitext(inp)[0]
            if mol_files:
                for mol in mol_files:
                    mol_no_suffix = os.path.splitext(mol)[0]
                    if other_files:
                        for other in other_files:
                            other_no_suffix = os.path.splitext(other)[0]
                            if inp_no_suffix == mol_no_suffix and mol_no_suffix == other_no_suffix:
                                output_no_suffix = '%s' % (inp_no_suffix,)
                            else:
                                output_no_suffix = '%s_%s_%s' % (inp_no_suffix, mol_no_suffix, other_no_suffix)
                            outputs_no_suffix.append(output_no_suffix)
                            messages.append('\nrunning test: %s %s %s\n' % (inp_no_suffix, mol_no_suffix, other_no_suffix))
                            commands.append(launcher + ' %s %s %s' % (inp_no_suffix, mol_no_suffix, other_no_suffix))
                    else:
                        if inp_no_suffix == mol_no_suffix:
                            output_no_suffix = '%s' % (inp_no_suffix,)
                        else:
                            output_no_suffix = '%s_%s' % (inp_no_suffix, mol_no_suffix)
                        outputs_no_suffix.append(output_no_suffix)
                        messages.append('\nrunning test: %s %s\n' % (inp_no_suffix, mol_no_suffix))
                        commands.append(launcher + ' %s %s' % (inp_no_suffix, mol_no_suffix))
            else:
                outputs_no_suffix.append('%s' % (inp_no_suffix,))
                messages.append('\nrunning test: %s\n' % (inp_no_suffix,))
                commands.append(launcher + ' %s' % (inp_no_suffix,))
        for output_no_suffix, message, command in zip(outputs_no_suffix, messages, commands):
            try:
                sys.stdout.write(message)
                runtest.TestRun.execute(self,
                                        command=command,
                                        stdout_file_name = '%s.stdout' % output_no_suffix,
                                        accepted_errors=accepted_errors)
                if f is None:
                    sys.stdout.write('finished (no reference)\n')
                else:
                    try:
                        # f is a suffix-filter dictionary
                        for suffix in f:
                            out = '%s.%s' % (output_no_suffix, suffix)
                            f[suffix].check(self.work_dir, '%s' % out, 'result/%s' % out, self.verbose)
                        sys.stdout.write('passed\n')
                    except IOError as e:
                        write_stderr(self.log, 'ERROR: could not open file %s\n' % e.filename)
                        sys.exit(-1)
                    except runtest.TestFailedError as e:
                        write_stderr(self.log, str(e))
                        self.return_code += 1
                    except runtest.BadFilterError as e:
                        write_stderr(self.log, str(e))
                        sys.exit(-1)
                    except runtest.FilterKeywordError as e:
                        write_stderr(self.log, str(e))
                        sys.exit(-1)
            except runtest.AcceptedError as e:
                sys.stdout.write(str(e))
            except runtest.SubprocessError as e:
                write_stderr(self.log, str(e))
                sys.exit(-1)
