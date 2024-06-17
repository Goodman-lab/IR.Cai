#!/usr/bin/python

import sys
import os
from csv import writer
from math import exp, sqrt
from pathlib import Path
import logging
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class IRAnalysis:
    def __init__(self, config_file):
        self.config_file = config_file
        self.molecule_title = Path(config_file).stem
        self.minimum_wavenumber = 1000.0
        self.maximum_wavenumber = 2000.0
        self.hwhm_lor_broad = 5.0
        self.temperature = 300.0
        self.R_in_hartrees = 0.0000031668
        self.defined_scaling_factor = 0.975
        self.max_scale_factor = 1.0
        self.min_scale_factor = 0.95
        self.optimise_scaling_factor = False
        self.extreme_sf_warning = 0.01
        self.max_energy_difference = 0.0001
        self.max_vcd_frequency_difference = 2.0
        self.boltzmann_analysis = True
        self.boltzmann_cutoff = 10.0
        self.print_graph = True
        self.reverse_xaxis = True
        self.print_spectra_csv = False
        self.print_scaling_factor_csv = False
        self.print_summary_csv = False
        self.calc_files = []
        self.ir_expt_file = ""
        self.expt_wavenums = []
        self.ir_expt_signals = []
        self.calc_energies = []
        self.dft = ""
        self.basis_set = ""
        self.optimise = ""
        self.calc_files_unique = []
        self.calc_energies_unique = []
        self.calc_boltzmann = []
        self.calc_boltzmann_unique = []
        self.calc_spectrum_wavenumbers_boltz = []
        self.calc_spectrum_wnintensity_boltz = []
        self.count_energy_cutoff_rejects = 0
        self.save_path = ""

        # Setup logger
        logging.basicConfig(
            level=logging.INFO,
            format='%(message)s',
            handlers=[
                logging.FileHandler(f"{self.molecule_title}.log"),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger()

    def parse_config(self):
        section_parsers = {
            "<settings>": self.parse_settings,
            "<experiments>": self.parse_experiments,
            "<calculations>": self.parse_calculations,
        }

        with open(self.config_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith("#"):
                    continue
                for section, parser in section_parsers.items():
                    if line.lower().startswith(section):
                        parser(f)
                        break

    def parse_settings(self, f):
        for line in f:
            line = line.strip()
            if line.lower().startswith("</settings>"):
                break
            self.parse_setting(line)

    def parse_setting(self, line):
        settings_map = {
            "title:": ("molecule_title", str),
            "broadening": ("hwhm_lor_broad", float),
            "minimum_wavenumber": ("minimum_wavenumber", float),
            "maximum_wavenumber": ("maximum_wavenumber", float),
            "temperature": ("temperature", float),
            "defined_scaling_factor": ("defined_scaling_factor", float),
            "max_scale_factor": ("max_scale_factor", float),
            "min_scale_factor": ("min_scale_factor", float),
            "optimise_scaling_factor": ("optimise_scaling_factor", lambda x: x.lower() != "false"),
            "save_path": ("save_path", str),
            "boltzmann_analysis": ("boltzmann_analysis", lambda x: x.lower() != "false"),
            "print_graph": ("print_graph", lambda x: x.lower() != "false"),
            "reverse_xaxis": ("reverse_xaxis", lambda x: x.lower() != "false"),
            "boltzmann_cutoff": ("boltzmann_cutoff", float),
            "extreme_sf_warning": ("extreme_sf_warning", float),
            "max_unique_energy_difference": ("max_energy_difference", float),
            "max_unique_peak_frequency_difference": ("max_vcd_frequency_difference", float),
        }

        # Check for print_csv separately since it can have multiple options
        if "print_csv" in line.lower():
            if "spectra" in line.lower():
                self.print_spectra_csv = True
            if "scaling_factor" in line.lower():
                self.print_scaling_factor_csv = True
            if "summary" in line.lower():
                self.print_summary_csv = True
            return

        for key, (attr, attr_type) in settings_map.items():
            if key in line.lower():
                setattr(self, attr, attr_type(line.split()[1]))
                return

    def parse_experiments(self, f):
        for line in f:
            line = line.strip()
            if line.lower().startswith("</experiments>"):
                break
            self.ir_expt_file = line

    def parse_calculations(self, f):
        for line in f:
            line = line.strip()
            if line.lower().startswith("</calculations>"):
                break
            self.calc_files.append(line)

    def get_ir_data(self, spec_file, wav_min, wav_max):
        lst_wnum = []
        lst_intensity = []
        with open(spec_file) as data_file:
            for line in data_file:
                linespace = line.replace(";", " ")
                try:
                    wnum = float(linespace.split()[0].strip())
                    intensity = float(linespace.split()[1].strip())
                    if wav_min <= wnum <= wav_max:
                        lst_wnum.append(wnum)
                        lst_intensity.append(intensity)
                except:
                    pass
        return lst_wnum, lst_intensity

    def get_calc_ir_data(self, calc_file, wav_min, wav_max, calc_file_suffix):
        calc_data_wavenumber = []
        calc_data_intensity = []
        try:
            with open(calc_file + calc_file_suffix) as f:
                for line in f:
                    if "Frequencies" in line:
                        for i in range(2, len(line.split())):
                            calc_data_wavenumber.append(float(line.split()[i].strip()))
                    if "IR Inten" in line:
                        for i in range(3, len(line.split())):
                            calc_data_intensity.append(float(line.split()[i].strip()))
        except IOError as e:
            self.logger.error(f"File open error for: {calc_file}{calc_file_suffix}")
        return calc_data_wavenumber, calc_data_intensity

    def get_calc_energy(self, calc_files, calc_file_suffix):
        calc_energy = []
        optimise = "Single point"
        dft = "?"
        basis_set = "?"
        for filename in calc_files:
            check_energy_calc = len(calc_energy)
            with open(filename + calc_file_suffix) as f:
                for line in f:
                    if "Sum of electronic and thermal Free Energies" in line:
                        calc_energy.append(float(line.split()[-1]))
                    elif "optimizer" in line:
                        optimise = "Optimized"
                    elif "freq=VCD" in line:
                        basis_set = line.split()[1].split("/")[1]
                        dft = line.split()[1].split("/")[0]
                if len(calc_energy) == check_energy_calc:
                    calc_energy.append(0.000)
        return calc_energy, dft, basis_set, optimise

    def lorentzian(self, expt_wavenums, calc_spectrum_wavenumbers, calc_spectrum_wnintensity, hwhm_lor_broad, scaling_factor):
        lor_data = []
        for wavenumber in expt_wavenums:
            lor_data.append(0.0)
            for i_peak in range(len(calc_spectrum_wavenumbers)):
                x_value = 2.0 * (wavenumber - calc_spectrum_wavenumbers[i_peak] * scaling_factor) / hwhm_lor_broad
                lor_data[-1] += calc_spectrum_wnintensity[i_peak] / (1 + x_value * x_value)
        return lor_data

    def match_score_calc(self, calc_signals, expt_signals):
        sum_alfacalc = 0.0
        sum_calc2 = 0.0
        sum_a2 = 0.0
        for i in range(len(calc_signals)):
            sum_alfacalc += calc_signals[i] * expt_signals[i]
            sum_calc2 += calc_signals[i] * calc_signals[i]
            sum_a2 += expt_signals[i] * expt_signals[i]
        return [sum_alfacalc / sqrt(sum_calc2 * sum_a2), sum_alfacalc]

    def print_graph_files(self, ir_cai, scale_factor):
        fig, ax1 = plt.subplots(1, 1)
        if self.reverse_xaxis:
            ax1.invert_xaxis()
        max_ir = max(self.ir_expt_signals)
        max_ir_calc = max(self.calc_signals)
        calc_scale_factor = max_ir_calc / max_ir
        scale_ir_calc_signals = [sig / calc_scale_factor for sig in self.calc_signals]
        ax1.set_title(f'IR.Cai: {ir_cai:5.3}; Scale Factor {scale_factor:6.4}\n{self.config_file.split(".")[0]}')
        ax1.plot(self.expt_wavenums, self.ir_expt_signals, label='IR', color='blue')
        ax1.plot(self.expt_wavenums, scale_ir_calc_signals, label='IR calc', color='cyan')
        ax1.legend()
        ax1.set_xlabel('wavenumbers')
        fig.savefig(self.config_file.split(".")[0] + '.pdf')
        plt.close(fig)

    def read_experimental_data(self):
        self.expt_wavenums, self.ir_expt_signals = self.get_ir_data(self.ir_expt_file, self.minimum_wavenumber, self.maximum_wavenumber)

    def read_calculation_files(self):
        self.update_calc_files_list()
        
        self.calc_energies, self.dft, self.basis_set, self.optimise = self.get_calc_energy(self.calc_files, '')
        self.sort_files_by_energy()

        filename = self.calc_files[0]
        calc_spectrum_wavenumbers, calc_spectrum_wnintensity = self.get_calc_ir_data(
            filename, self.minimum_wavenumber, self.maximum_wavenumber, ''
        )
        self.calc_spectrum_wavenumbers_boltz = calc_spectrum_wavenumbers
        self.calc_spectrum_wnintensity_boltz = calc_spectrum_wnintensity
        number_of_signals = [len(calc_spectrum_wavenumbers)]

        if len(self.calc_files) > 1 and self.boltzmann_analysis:
            self.perform_boltzmann_analysis(number_of_signals)

        self.log_summary()

    def update_calc_files_list(self):
        if len(self.calc_files) == 1 and Path(self.calc_files[0]).is_dir():
            path_list = list(Path(self.calc_files[0]).glob('*.log'))
            self.calc_files = [os.path.splitext(str(filename))[0] + ".log" for filename in path_list]
            if not self.calc_files:
                self.logger.error("Error: no files")
                exit()
            self.logger.info("Number of files: %d", len(self.calc_files))

    def sort_files_by_energy(self):
        pairs = list(zip(self.calc_energies, self.calc_files))
        sort_pairs = sorted(pairs)
        self.calc_energies = [p[0] for p in sort_pairs]
        self.calc_files = [p[1] for p in sort_pairs]

    def perform_boltzmann_analysis(self, number_of_signals):
        self.calc_spectrum_wavenumbers_boltz = []
        self.calc_spectrum_wnintensity_boltz = []

        for i in range(len(self.calc_files)):
            boltzmann_factor = exp(-(self.calc_energies[i] - self.calc_energies[0]) / self.R_in_hartrees / self.temperature)
            self.calc_boltzmann.append(boltzmann_factor)
            
            filename = self.calc_files[i]
            calc_spectrum_wavenumbers, calc_spectrum_wnintensity = self.get_calc_ir_data(
                filename, self.minimum_wavenumber, self.maximum_wavenumber, ''
            )
            if self.is_within_energy_cutoff(i):
                self.update_signals_and_structures(i, calc_spectrum_wavenumbers, calc_spectrum_wnintensity, number_of_signals)

    def is_within_energy_cutoff(self, i):
        return self.calc_energies[i] - self.calc_energies[0] < self.boltzmann_cutoff / 2625.5

    def update_signals_and_structures(self, i, calc_spectrum_wavenumbers, calc_spectrum_wnintensity, number_of_signals):
        if len(calc_spectrum_wavenumbers) > 0:
            number_of_signals.append(len(calc_spectrum_wavenumbers))
            if self.is_unique_structure(i, calc_spectrum_wavenumbers):
                self.store_unique_structure(i, calc_spectrum_wavenumbers, calc_spectrum_wnintensity)

    def is_unique_structure(self, i, calc_spectrum_wavenumbers):
        for j in range(i):
            if abs(self.calc_energies[j] - self.calc_energies[i]) < self.max_energy_difference:
                if self.compare_spectra(j, calc_spectrum_wavenumbers):
                    return False
        return True

    def compare_spectra(self, j, calc_spectrum_wavenumbers):
        test_filename = self.calc_files[j]
        test_calc_spectrum_wavenumbers, _ = self.get_calc_ir_data(
            test_filename, self.minimum_wavenumber, self.maximum_wavenumber, ''
        )
        if len(test_calc_spectrum_wavenumbers) != len(calc_spectrum_wavenumbers):
            return False
        for k in range(len(calc_spectrum_wavenumbers)):
            if abs(test_calc_spectrum_wavenumbers[k] - calc_spectrum_wavenumbers[k]) > self.max_vcd_frequency_difference:
                return False
        return True

    def store_unique_structure(self, i, calc_spectrum_wavenumbers, calc_spectrum_wnintensity):
        boltzmann_factor = exp(-(self.calc_energies[i] - self.calc_energies[0]) / self.R_in_hartrees / self.temperature)
        self.calc_energies_unique.append(self.calc_energies[i])
        self.calc_files_unique.append(self.calc_files[i])
        self.calc_boltzmann_unique.append(boltzmann_factor)
        for j in range(len(calc_spectrum_wavenumbers)):
            self.calc_spectrum_wavenumbers_boltz.append(calc_spectrum_wavenumbers[j])
            self.calc_spectrum_wnintensity_boltz.append(calc_spectrum_wnintensity[j] * boltzmann_factor)

    def log_summary(self):
        self.logger.info("%d files rejected by energy cutoff; %d duplicate files removed", 
                        self.count_energy_cutoff_rejects, 
                        len(self.calc_energies) - len(self.calc_energies_unique) - self.count_energy_cutoff_rejects)

        self.logger.info("Unique Calculated Structures")
        for i in range(len(self.calc_energies_unique)):
            self.logger.info("     Energy: %10.6f hartrees, %6.3f kJ/mol, Boltzmann Factor: %5.3f %s", 
                            self.calc_energies_unique[i], 
                            (self.calc_energies_unique[i] - self.calc_energies_unique[0]) * 2625.8, 
                            self.calc_boltzmann_unique[i], 
                            self.calc_files_unique[i])

        if self.boltzmann_analysis:
            self.logger.info("Using all %d unique conformations within energy cut-off", len(self.calc_files_unique))
        else:
            self.logger.info("Using only lowest energy conformation in analysis: %s Energy: %f", 
                            self.calc_files_unique[0], 
                            self.calc_energies_unique[0])
        self.logger.info("")

    def calculate_signals(self):
        self.calc_signals = self.lorentzian(self.expt_wavenums, self.calc_spectrum_wavenumbers_boltz, self.calc_spectrum_wnintensity_boltz, self.hwhm_lor_broad, self.defined_scaling_factor)

    def perform_scaling_factor_analysis(self):
        best_sf_ir_cai = 0.0
        best_sf_ir = 1.0

        if self.optimise_scaling_factor:
            scale_factor_range = self.max_scale_factor - self.min_scale_factor
            scale_scale = 20.0 / scale_factor_range if scale_factor_range > 0.001 else 20000

            for scaling_factor_scale in range(int(self.min_scale_factor * scale_scale), int(self.max_scale_factor * scale_scale + 1.0)):
                scaling_factor = float(scaling_factor_scale) / scale_scale
                if scaling_factor < self.min_scale_factor:
                    scaling_factor = self.min_scale_factor
                if scaling_factor > self.max_scale_factor:
                    scaling_factor = self.max_scale_factor

                self.calc_signals = self.lorentzian(self.expt_wavenums, self.calc_spectrum_wavenumbers_boltz, self.calc_spectrum_wnintensity_boltz, self.hwhm_lor_broad, scaling_factor)
                ir_cai = self.match_score_calc(self.calc_signals, self.ir_expt_signals)[0]

                if best_sf_ir_cai < ir_cai:
                    best_sf_ir_cai = ir_cai
                    best_sf_ir = scaling_factor

        # Calculate match score using the defined scaling factor
        self.calc_signals = self.lorentzian(self.expt_wavenums, self.calc_spectrum_wavenumbers_boltz, self.calc_spectrum_wnintensity_boltz, self.hwhm_lor_broad, self.defined_scaling_factor)
        defined_sf_ir_cai = self.match_score_calc(self.calc_signals, self.ir_expt_signals)[0]

        output_file = open(self.config_file.split(".")[0] + "_output.csv", 'w')
        print(f"wavenumbers,ir_data,ir_calc,SF {self.defined_scaling_factor:5.3f}", file=output_file)
        for i in range(len(self.expt_wavenums)):
            print(f"{self.expt_wavenums[i]}, {self.ir_expt_signals[i]}, {self.calc_signals[i]}", file=output_file)
        output_file.close()

        return best_sf_ir, best_sf_ir_cai, defined_sf_ir_cai

    def log_results(self, best_sf_ir, best_sf_ir_cai, defined_sf_ir_cai):
        self.logger.info(f'Best SF:      SF = {best_sf_ir:.4f}, IR.Cai = {best_sf_ir_cai:.4f}')

        self.logger.info(f'Defined SF:   SF = {self.defined_scaling_factor:.4f}, IR.Cai = {defined_sf_ir_cai:.4f}')

        if self.save_path:
            with open(self.save_path, 'a') as f:
                writer_object = writer(f)
                cwd = os.getcwd()
                molecule = os.path.basename(cwd)
                result = [molecule, defined_sf_ir_cai]
                writer_object.writerow(result)

    def analyse(self):
        self.parse_config()
        output_string = (
            "##########################################\n"
            "##          IR.Cai analysis             ##\n"
            "##    University of Cambridge, 2024     ##\n"
            "##########################################\n")
        self.logger.info(output_string)
        self.logger.info(self.config_file + "\n")

        self.read_experimental_data()
        self.read_calculation_files()
        self.calculate_signals()
        best_sf_ir, best_sf_ir_cai, defined_sf_ir_cai = self.perform_scaling_factor_analysis()
        self.log_results(best_sf_ir, best_sf_ir_cai, defined_sf_ir_cai)

        if self.print_graph:
            self.print_graph_files(defined_sf_ir_cai, self.defined_scaling_factor)

if __name__ == '__main__':
    config_file = sys.argv[1]
    analysis = IRAnalysis(config_file)
    analysis.analyse()