#!/usr/bin/env python3
# -*-coding:utf-8-*-

__version__ = "v2.0.2"
__author__ = "Jia Liu"

__str__ = """
********** HELP MANUSCRIPT **********

[INTRODUCTION]
  Calculate effective masses of holes and electrons automatically. (for DMol3 and CASTEP)
  Note that "Band structure" calculation and analyzation should be performed before running this script.
  This script is released under GPL v3.0 license.
  Current version: %s
[AUTHOR] %s
[WEBSITE] https://github.com/liujiacode/EffectiveMass

[USAGE]
  >> python EffectiveMass.py -p "project_directory" [-s "md_xsd_file"] [-c "bs_xcd_file"] [-h] [-v]

[OUTPUT]
  This script would create a result file named "XXX Effective Mass.out" in project directory.

[OPTIONS]
-p (or --project):
  This is a required option to indicate the project directory.
  The name of project should follow the naming scheme.
  For example, >> python EffectiveMass.py -p "/path/to/XXX DMol3 GeomOpt"
  See also "project_directory" parameter.

-s (or --xsd):
  This is an optional option to indicate the model file with xsd extending.
  Default loading "XXX.xsd" in project directory if model file is not specified.
  For example, >> python EffectiveMass.py -p "/path/to/XXX DMol3 GeomOpt" -s "XXX(2).xsd"
  See also "md_xsd_file" parameter.

-c (or --xcd):
  This is an optional option to indicate the model file with xsd extending.
  Default loading "XXX Band Structure.xcd" in project directory if band structure file is not specified.
  For example, >> python EffectiveMass.py -p "/path/to/XXX DMol3 GeomOpt" -c "XXX Band Structure(2).xcd"
  See also "bs_xcd_file" parameter.

-h (or --help):
  This is an optional option without parameters to view the help manuscript.
  For example, >> python EffectiveMass.py -h

-v (or --version):
  This is an optional option without parameters to check the version of EffectiveMass.py script.
  For example, >> python EffectiveMass.py -v

[PARAMETERS]
project directory ("project_directory"):
  The full path of project directory.
  This script supports project calculated by DMol3 (or CASTEP) code with GeomOpt (or Energy) method.
  So the project name should contains both code and method name as the ending words.
  For example, "XXX DMol3 GeomOpt" is accepted, but "XXX", "DMol3 GeomOpt", "XXX DMol3" or "XXX GeomOpt", etc. are illegal.
  The project directory should be indicated manually by using "-p" option in this script.
  See also "-p (or --project)" option.

model file with xsd extending ("md_xsd_file"):
  The full name of model file.
  The model file should be contained in the project directory and the full name is default as "XXX.xsd".
  This file is generated automatically by running DMol3 (or CASTEP) Calculation.
  If you want to specify another model file such as "XXX(2).xsd", please use "-s" option.
  See also "-s (or --xsd)" option.

band structure file with xcd extending ("bs_xcd_file"):
  The full name of band structure file.
  The band structure file should be contained in the project directory and the full name is default as "XXX Band Structure.xcd".
  This file is generated automatically by running DMol3 (or CASTEP) Analysis.
  If you want to specify another band structure file such as "XXX Band Structure(2).xcd", please use "-c" option.
  See also "-c (or --xcd)" option.

**************** END ****************
""" % (__version__, __author__)

"""
v2.0.0 UPDATE:
  The calculations of hole and electron effective masses are confused in v1.x.x, 
  and they have been corrected in v2.0.0.

v2.0.1 UPDATE:
  Sometimes x index in band structure end up with 0.99999... instead of 1.0, 
  which would lead to error "Unknown length of x axis" in v2.0.0 or older.
"""

import math
import os
import sys
from lxml import html
import re
import getopt


class Vector(object):
    def __init__(self, vector):
        if isinstance(vector, Vector):
            self._x, self._y, self._z = vector.values()
        elif isinstance(vector, tuple):
            if len(vector) == 3:
                self._x, self._y, self._z = (float(i) for i in vector)
            else:
                raise ValueError("Invalid input tuple: %s." % vector)
        else:
            raise ValueError("Invalid input vector: %s." % vector)

    def values(self):
        return self._x, self._y, self._z

    def __repr__(self):
        return str(self.values())

    def __getitem__(self, item):
        if isinstance(item, int) and 0 <= item <= 2:
            return self.values()[item]
        elif isinstance(item, str) and item in ('x', 'y', 'z'):
            return {'x': self._x, 'y': self._y, 'z': self._z}[item]
        else:
            raise ValueError("Invalid item: {}.".format(item))

    def __add__(self, other):
        if isinstance(other, Vector):
            return Vector(tuple(map(lambda x, y: x + y, self.values(), other.values())))
        else:
            raise ValueError("Invalid vector for add: {}.".format(other))

    def __sub__(self, other):
        if isinstance(other, Vector):
            return Vector(tuple(map(lambda x, y: x - y, self.values(), other.values())))
        else:
            raise ValueError("Invalid vector for sub: {}.".format(other))

    def __mul__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return Vector(tuple(map(lambda x: x * other, self.values())))
        else:
            raise ValueError("Invalid number for mul: {}.".format(other))

    def __truediv__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            if other == 0.0:
                raise ValueError("Invalid truediv for zero.")
            else:
                return Vector(tuple(map(lambda x: x / other, self.values())))
        else:
            raise ValueError("Invalid number for truediv: {}.".format(other))

    def __floordiv__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            if other == 0.0:
                raise ValueError("Invalid truediv for zero.")
            else:
                return Vector(tuple(map(lambda x: x // other, self.values())))
        else:
            raise ValueError("Invalid number for floordiv: {}.".format(other))

    def number_product(self, other):
        if isinstance(other, Vector):
            return sum(tuple(map(lambda x, y: x * y, self.values(), other.values())))
        else:
            raise ValueError("Invalid vector for number product: {}.".format(other))

    def vector_product(self, other):
        if isinstance(other, Vector):
            if self == Vector((0, 0, 0)) or other == Vector((0, 0, 0)):
                raise ValueError("Invalid vector product for zero vector.")
            else:
                return Vector((self._y * other['z'] - other['y'] * self._z,
                               self._z * other['x'] - other['z'] * self._x,
                               self._x * other['y'] - other['x'] * self._y))
        else:
            raise ValueError("Invalid vector for vector product: {}.".format(other))

    def __eq__(self, other):
        if isinstance(other, Vector):
            if self.values() == other.values():
                return True
            else:
                return False
        else:
            raise ValueError("Invalid vector for abs_equal: {}.".format(other))

    def __ne__(self, other):
        if isinstance(other, Vector):
            if self.values() == other.values():
                return False
            else:
                return True
        else:
            raise ValueError("Invalid vector for ne: {}.".format(other))

    def length(self):
        return math.sqrt(sum(map(lambda x: x ** 2, self.values())))


class EffectiveMass(object):
    def __init__(self, project):
        # calculating precision
        self._cal_precision = 1e-3
        # Processing precision
        self._pro_precision = 1e-6

        if os.path.isdir(project):
            self._project = project
        else:
            raise ValueError("Project is not exist: %s." % project)

        self._project_name = os.path.basename(self._project)
        if len(self._project_name.split(' ')) < 3:
            raise ValueError("Can not identify the calculation code and method: %s." % self._project_name)
        self._name = ' '.join(self._project_name.split(' ')[:-2])
        self._code = self._project_name.split(' ')[-2]
        self._method = self._project_name.split(' ')[-1]
        if self._code not in ["DMol3", "CASTEP"]:
            raise ValueError("Calculation code is not supported: %s." % self._code)
        if self._method not in ["GeomOpt", "Energy"]:
            raise ValueError("Calculation method is not supported: %s." % self._project_name)

        output_file = self._project + os.sep + self._name + " Effective Mass"
        file_index = 1
        if os.path.isfile(output_file + ".out"):
            while os.path.isfile(output_file + '(%s).out' % str(file_index)):
                file_index += 1
            output_file = output_file + '(%s)' % str(file_index)
        self._output_file = output_file + ".out"
        with open(self._output_file, 'w') as o_file:
            o_file.write("Electron and hole effective mass calculated by %s with argv %s.\n" % (sys.path[0] + os.sep + sys.argv[0], sys.argv[1:]))
            o_file.write("EffectiveMass.py version: %s.\n" % __version__)
            o_file.write("Project path: %s.\n" % os.path.abspath(self._project))
            o_file.write("Calculation code: %s\n" % self._code)
            o_file.write("Calculation method: %s\n" % self._method)
            o_file.write("\nPlease check the calculation results carefully!\n")

    def output(self, string):
        with open(self._output_file, 'a') as o_file:
            o_file.write(str(string) + '\n')

    def load(self, xsd_file="default", bs_xcd_file="default"):
        if xsd_file == "default":
            xsd_file = self._project + os.sep + self._name + ".xsd"
        else:
            xsd_file = self._project + os.sep + str(xsd_file)
        if os.path.isfile(xsd_file):
            xsd_data = html.fromstring(open(xsd_file, 'r').read().encode("utf-8"))
        else:
            raise ValueError("Model file is not exist: %s." % xsd_file)

        if bs_xcd_file == "default":
            bs_xcd_file = self._project + os.sep + self._name + " Band Structure.xcd"
        else:
            bs_xcd_file = self._project + os.sep + str(bs_xcd_file)
        if os.path.isfile(bs_xcd_file):
            bs_xcd_data = html.fromstring(open(bs_xcd_file, 'r').read().encode("utf-8"))
        else:
            raise ValueError("Band structure file is not exist: %s." % bs_xcd_file)

        # Loading real lattice.
        space_group = xsd_data.xpath("//spacegroup")
        if not len(space_group):
            raise ValueError("Space group is not found in model file: %s." % xsd_file)
        if len(space_group) > 1:
            raise ValueError("Too many space groups record in model file: %s" % xsd_file)
        a_vector = Vector(tuple(float(i) for i in space_group[0].xpath("./@avector")[0].split(',')))
        b_vector = Vector(tuple(float(i) for i in space_group[0].xpath("./@bvector")[0].split(',')))
        c_vector = Vector(tuple(float(i) for i in space_group[0].xpath("./@cvector")[0].split(',')))
        cell_vectors = (a_vector, b_vector, c_vector)

        # Loading k points and k lines from xsd_file.
        k_path = xsd_data.xpath("//kpath")
        if not len(k_path):
            raise ValueError("K path is not found in model file: %s." % xsd_file)
        if len(k_path) > 1:
            raise ValueError("Too many K paths record in model file: %s." % xsd_file)
        k_path_objects = list(int(i) for i in k_path[0].xpath("./@objects")[0].split(','))

        k_points = {}
        k_names = {}
        for k_point in xsd_data.xpath("//kpoint"):
            k_point_id = int(k_point.xpath("./@id")[0])
            k_point_name = k_point.xpath("./@name")[0]
            k_point_HKL = Vector(tuple(float(i) for i in k_point.xpath("./@hkl")[0].split(',')))
            k_points[k_point_name] = k_point_HKL
            k_names[k_point_id] = k_point_name
        if not k_points:
            raise ValueError("K points are not found in model file: %s." % xsd_file)

        k_lines = []
        for object_id in k_path_objects:
            k_line = xsd_data.xpath("//kline[@id=%s]" % object_id)
            if len(k_line) == 1:
                k_line_connects = tuple(int(i) for i in k_line[0].xpath("./@connects")[0].split(','))
                k_lines.append(tuple(k_names[i] for i in k_line_connects))
            else:
                raise ValueError("Can not identify K line with k path id: %s" % object_id)

        k_lines = tuple(k_lines)
        if not k_lines:
            raise ValueError("K lines are not found in model file: %s." % xsd_file)

        # Loading band structures.
        series_2d = bs_xcd_data.xpath("//series_2d")
        if not len(series_2d):
            raise ValueError("Series_2D are not found in band structure file: %s." % bs_xcd_file)
        if len(series_2d) > 1:
            raise ValueError("Too many Series_2D record in band structure file: %s." % bs_xcd_file)

        numpoints = int(series_2d[0].xpath("@numpoints")[0])
        point_2d = series_2d[0].xpath("./point_2d")
        if len(point_2d) != numpoints:
            raise ValueError("Unknown total number of point_2D: %s." % str(len(point_2d)))

        # Valence band and conduction band data share one x axis.
        valence_band = tuple()
        conduction_band = tuple()
        x_axis = []
        next_conduction_band = False
        vb_index = 0
        cb_index = 0
        for each_point in point_2d:
            x, y = [float(i) for i in each_point.xpath("@xy")[0].split(',')]
            x_axis.append(x)
            if abs(x - 1) < self._pro_precision:
                break
        x_axis_len = len(x_axis)
        # Each group of band data has one waste record like "<POINT_2D XY="0,1e+308"/>".
        # So it has to add 1 when splitting band data into groups.
        rx_axis_len = x_axis_len + 1
        if numpoints % rx_axis_len != 0:
            raise ValueError("Unknown length of x axis: %s." % str(rx_axis_len))

        for band_index in range(int(numpoints / rx_axis_len)):
            y_axis = []
            for point_index in range(x_axis_len):
                x, y = [float(i) for i in
                        point_2d[rx_axis_len * band_index + point_index].xpath("@xy")[0].split(',')]
                if x == x_axis[point_index]:
                    y_axis.append(y)
                else:
                    i_x_index = rx_axis_len * band_index + point_index
                    raise ValueError("Invalid x of point: %s." % (str(i_x_index)))

            dis_x, dis_y = [float(i) for i in point_2d[rx_axis_len * band_index + x_axis_len].xpath("@xy")[0].split(',')]
            if dis_x != 0:
                i_dis_x = rx_axis_len * band_index + x_axis_len
                raise ValueError("Invalid dis_x of last discarded point: %s." % str(i_dis_x))

            if next_conduction_band:
                cb_index = band_index
                conduction_band = tuple(y_axis)
                next_conduction_band = False
            if abs(max(y_axis)) < 0.0001:
                vb_index = band_index
                valence_band = tuple(y_axis)
                next_conduction_band = True

        self.output('')
        self.output("++++++++++++ Loading ++++++++++++")
        self.output('Find valence band at index: %s' % vb_index)
        self.output('Find conduction band at index: %s' % cb_index)

        x_axis = tuple(x_axis)
        if not valence_band:
            raise ValueError("Loading valence band data failed.")
        if not conduction_band:
            raise ValueError("Loading conduction band data failed.")

        # Loading other paras.
        cart_2d_plot = bs_xcd_data.xpath("//cartesian_2d_plot")
        if len(cart_2d_plot) != 1:
            raise ValueError("Too many(little) 'CARTESIAN_2D_PLOT' patterns: %s." % bs_xcd_file)
        plot_title = cart_2d_plot[0].xpath("./title/@text")[0]
        title_re = re.compile(r"^%s Band Structure\\n Band gap is (\d+\.\d+) (Ha|eV)$" % self._code)
        title_match = title_re.match(plot_title)
        if title_match:
            band_gap_num, band_gap_unit = title_match.groups()
            band_gap_num = float(band_gap_num)
        else:
            raise ValueError("Title is mismatching: %s." % plot_title)
        calculated_band_gap_num = min(conduction_band) - max(valence_band)
        if abs(calculated_band_gap_num - band_gap_num) > self._cal_precision:
            raise ValueError("Invalid calculated band gap value: %s." % calculated_band_gap_num)

        # Loading k path labels (high symmetry points in band structure files).
        num_labels = int(cart_2d_plot[0].xpath("./cartesian_axes_2d/axis_tick_labels/@numlabels")[0])
        axis_tick_labels = cart_2d_plot[0].xpath("./cartesian_axes_2d/axis_tick_labels/axis_tick_label")
        if len(axis_tick_labels) != num_labels:
            raise ValueError("Unknown axis tick labels: %s." % num_labels)
        k_path_labels = []
        for axis_tick_label in axis_tick_labels:
            position = float(axis_tick_label.xpath("./@position")[0])
            text = axis_tick_label.xpath("./@text")[0]
            if text in k_points:
                k_path_labels.append((text, position))
            else:
                raise ValueError("Point %s: %s is not found in k points list." % (text, position))

        k_path_labels = tuple(k_path_labels)
        if not k_path_labels:
            raise ValueError("Loading k path labels data failed.")

        # The sequence of k_lines should match with k_path_labels.
        for index in range(len(k_lines)):
            start = k_lines[index][0]
            end = k_lines[index][1]
            if start != k_path_labels[index][0] or end != k_path_labels[index + 1][0]:
                raise ValueError("K lines and K path labels are mismatching.")

        k_path_paras = (k_points, k_lines, k_path_labels)
        band_paras = (band_gap_num, band_gap_unit, x_axis, valence_band, conduction_band)
        return cell_vectors, k_path_paras, band_paras

    def calculate(self, cell_vectors, k_path_paras, band_paras):
        k_points, k_lines, k_path_labels = k_path_paras
        band_gap_num, band_gap_unit, x_axis, valence_band, conduction_band = band_paras
        bohr_constant = 0.52917721092
        ha_constant = 27.2114
        # Calculate reciprocal vectors with bohr unit.
        bohr_cell_vectors = tuple(x / bohr_constant for x in cell_vectors)
        volume = (bohr_cell_vectors[0].vector_product(bohr_cell_vectors[1]).number_product(bohr_cell_vectors[2]))
        bohr_reciprocal_vectors = (
            bohr_cell_vectors[1].vector_product(bohr_cell_vectors[2]) * (2 * math.pi / volume),
            bohr_cell_vectors[2].vector_product(bohr_cell_vectors[0]) * (2 * math.pi / volume),
            bohr_cell_vectors[0].vector_product(bohr_cell_vectors[1]) * (2 * math.pi / volume),
        )

        # Calculate abs k points under reciprocal lattice (bohr_reciprocal_vectors).
        cart_k_points = {}
        for index in k_points.keys():
            cart_vector = Vector((0, 0, 0))
            for i in range(3):
                cart_vector = cart_vector + bohr_reciprocal_vectors[i] * k_points[index][i]
            cart_k_points[index] = cart_vector

        # Calculate the abs length of each k line.
        k_line_lengths = []
        for k_line in k_lines:
            start = cart_k_points[k_line[0]]
            end = cart_k_points[k_line[1]]
            k_line_lengths.append(Vector(end - start).length())

        # Calculate the abs labels of each k point along k lines.
        cart_k_path_labels = [k_path_labels[0]]
        for index in range(len(k_line_lengths)):
            abs_length = sum(k_line_lengths[:index + 1])
            cart_k_path_labels.append((k_path_labels[index + 1][0], abs_length))

        # Get the relative position of each k point in axis.
        kp_in_x_axis = []
        for k_path_label in k_path_labels:
            label = k_path_label[1]
            position = -1
            previous = position
            for x_index in range(len(x_axis)):
                if abs(label - x_axis[x_index]) < 0.0001:
                    position = x_index
                    break
            if position >= 0 and position > previous:
                kp_in_x_axis.append((k_path_label[0], position))
            else:
                raise ValueError("Can not find the position of K path label: %s." % str(k_path_label))
        if len(kp_in_x_axis) != len(k_path_labels):
            raise ValueError("Unknown K points in x axis: %s." % str(kp_in_x_axis))

        # Calculate the absolute x axis values.
        cart_x_axis = [-1] * len(x_axis)
        for index in range(len(kp_in_x_axis) - 1):
            start_label = kp_in_x_axis[index][0]
            end_label = kp_in_x_axis[index + 1][0]
            start_index = kp_in_x_axis[index][1]
            end_index = kp_in_x_axis[index + 1][1]
            start_value = x_axis[start_index]
            end_value = x_axis[end_index]
            if start_label == cart_k_path_labels[index][0]:
                cart_start_value = cart_k_path_labels[index][1]
            else:
                raise ValueError("The start label is mismatching: %s." % start_label)
            if end_label == cart_k_path_labels[index + 1][0]:
                cart_end_value = cart_k_path_labels[index + 1][1]
            else:
                raise ValueError("The end label is mismatching: %s." % end_label)
            rate = (cart_end_value - cart_start_value) / (end_value - start_value)
            for x_index in range(end_index - start_index + 1):
                abs_x_index = start_index + x_index
                abs_x_value = cart_start_value + (x_axis[abs_x_index] - start_value) * rate
                if cart_x_axis[abs_x_index] == -1:
                    cart_x_axis[abs_x_index] = abs_x_value
                else:
                    if abs(cart_x_axis[abs_x_index] - abs_x_value) > 0.0001:
                        raise ValueError("Cartesian x axis value is mismatching: %s." % str(cart_x_axis[abs_x_index]))

        for value in cart_x_axis:
            if value < 0:
                raise ValueError("Cartesian x axis value is less than zero: %s." % str(value))

        for index in range(len(cart_x_axis) - 1):
            if cart_x_axis[index] >= cart_x_axis[index + 1]:
                raise ValueError("Unknown Cartesian x axis value: %s." % str(cart_x_axis[index]))

        # This script can only calculate simple band structures.
        if kp_in_x_axis[0][0] != kp_in_x_axis[-1][0]:
            raise ValueError("The start K point label of band structure is not coincident with the end label.")

        pos_list = []
        for pos in kp_in_x_axis[:-1]:
            if pos[0] not in pos_list:
                pos_list.append(pos[0])
            else:
                raise ValueError("Too many K point in band structure: %s." % pos[0])

        # Effective mass of hole on valence band.
        hole_value = max(valence_band)
        hole_position = valence_band.index(hole_value)
        if hole_value != 0.0:
            raise ValueError("Maximum value of valence band is not zero.")

        # Two directions if hole_position is a k point.
        if hole_position in tuple(i[1] for i in kp_in_x_axis):
            hole_index = tuple(i[1] for i in kp_in_x_axis).index(hole_position)
            # First direction.
            if hole_index == 0:
                direction_1 = (kp_in_x_axis[-2], kp_in_x_axis[-1], len(cart_x_axis) - 1)
            else:
                direction_1 = (kp_in_x_axis[hole_index - 1], kp_in_x_axis[hole_index], hole_position)
            # Second direction.
            if hole_index == len(kp_in_x_axis) - 1:
                direction_2 = (kp_in_x_axis[0], kp_in_x_axis[1], 0)
            else:
                direction_2 = (kp_in_x_axis[hole_index], kp_in_x_axis[hole_index + 1], hole_position)
            hole_directions = (direction_1, direction_2)
        else:
            hole_directions = tuple()
            for hole_index in range(len(kp_in_x_axis) - 1):
                if kp_in_x_axis[hole_index][1] < hole_position < kp_in_x_axis[hole_index + 1][1]:
                    hole_directions = ((kp_in_x_axis[hole_index], kp_in_x_axis[hole_index + 1], hole_position),)
                    break
            if not hole_directions:
                raise ValueError("Hole direction is not found.")

        hole_effective_masses = []
        for hole_direction in hole_directions:
            if hole_direction[1][1] - hole_direction[0][1] < 4:
                raise ValueError("Number of points between hole direction: %s is less than four." % str(hole_direction))
            dir_name = "%s%s" % (hole_direction[0][0], hole_direction[1][0])
            start_pos = hole_direction[0][1]
            end_pos = hole_direction[1][1]
            point_pos = hole_direction[2]
            if point_pos == start_pos:
                p_1 = start_pos
                p_2 = start_pos + 1
                p_3 = start_pos + 2
            elif point_pos == end_pos:
                p_1 = end_pos - 2
                p_2 = end_pos - 1
                p_3 = end_pos
            else:
                p_1 = point_pos - 1
                p_2 = point_pos
                p_3 = point_pos + 1
            x_1 = cart_x_axis[p_1]
            x_2 = cart_x_axis[p_2]
            x_3 = cart_x_axis[p_3]
            y_1 = valence_band[p_1]
            y_2 = valence_band[p_2]
            y_3 = valence_band[p_3]
            valence_band_unit = band_gap_unit
            if valence_band_unit == 'eV':
                y_1 = y_1 / ha_constant
                y_2 = y_2 / ha_constant
                y_3 = y_3 / ha_constant
                valence_band_unit = 'Ha'
            if valence_band_unit != 'Ha':
                raise ValueError("Invalid valence band data unit: %s." % valence_band_unit)

            d_2 = (y_3 - y_2) / (x_3 - x_2)
            d_1 = (y_2 - y_1) / (x_2 - x_1)
            dd = (d_2 - d_1) / ((x_3 - x_1) / 2)
            mass = 1 / dd
            hole_effective_masses.append((dir_name, mass))

        if not hole_effective_masses:
            raise ValueError("Calculating hole effective mass failed, no hole effective mass is calculated.")

        self.output('')
        self.output("++++++++++++ Results ++++++++++++")
        self.output("m0 is the effective mass of a free electron.")
        self.output('')
        self.output('(1) Holes on valence band.')
        self.output('Max value on valence band: %s %s' % (hole_value, band_gap_unit))
        self.output("VBM position: %s" % hole_position)
        dir_content = ""
        for hole_dir in hole_directions:
            dir_content += "  %s(%s) -> %s(%s)\n" % (hole_dir[0][0], hole_dir[0][1], hole_dir[1][0], hole_dir[1][1])
        self.output("Hole direction(s):\n%s" % dir_content[:-1])
        mass_content = ""
        for hole_mass in hole_effective_masses:
            mass_content += "  Direction: %s, Mass: %s m0\n" % (hole_mass[0], hole_mass[1])
        self.output("Hole effective mass(es):\n%s" % mass_content[:-1])

        # Effective mass of electron on conduction band.
        elec_value = min(conduction_band)
        elec_position = conduction_band.index(elec_value)

        # Two directions if elec_position is a k point.
        if elec_position in tuple(i[1] for i in kp_in_x_axis):
            elec_index = tuple(i[1] for i in kp_in_x_axis).index(elec_position)
            # First direction.
            if elec_index == 0:
                direction_1 = (kp_in_x_axis[-2], kp_in_x_axis[-1], len(cart_x_axis) - 1)
            else:
                direction_1 = (kp_in_x_axis[elec_index - 1], kp_in_x_axis[elec_index], elec_position)
            # Second direction.
            if elec_index == len(kp_in_x_axis) - 1:
                direction_2 = (kp_in_x_axis[0], kp_in_x_axis[1], 0)
            else:
                direction_2 = (kp_in_x_axis[elec_index], kp_in_x_axis[elec_index + 1], elec_position)
            elec_directions = (direction_1, direction_2)
        else:
            elec_directions = tuple()
            for elec_index in range(len(kp_in_x_axis) - 1):
                if kp_in_x_axis[elec_index][1] < elec_position < kp_in_x_axis[elec_index + 1][1]:
                    elec_directions = ((kp_in_x_axis[elec_index], kp_in_x_axis[elec_index + 1], elec_position),)
                    break
            if not elec_directions:
                raise ValueError("Electron direction is not found.")

        elec_effective_masses = []
        for elec_direction in elec_directions:
            if elec_direction[1][1] - elec_direction[0][1] < 4:
                raise ValueError("Number of points between electron direction: %s is less than four." % str(elec_direction))
            dir_name = "%s%s" % (elec_direction[0][0], elec_direction[1][0])
            start_pos = elec_direction[0][1]
            end_pos = elec_direction[1][1]
            point_pos = elec_direction[2]
            if point_pos == start_pos:
                p_1 = start_pos
                p_2 = start_pos + 1
                p_3 = start_pos + 2
            elif point_pos == end_pos:
                p_1 = end_pos - 2
                p_2 = end_pos - 1
                p_3 = end_pos
            else:
                p_1 = point_pos - 1
                p_2 = point_pos
                p_3 = point_pos + 1
            x_1 = cart_x_axis[p_1]
            x_2 = cart_x_axis[p_2]
            x_3 = cart_x_axis[p_3]
            y_1 = conduction_band[p_1]
            y_2 = conduction_band[p_2]
            y_3 = conduction_band[p_3]
            conduction_band_unit = band_gap_unit
            if conduction_band_unit == 'eV':
                y_1 = y_1 / ha_constant
                y_2 = y_2 / ha_constant
                y_3 = y_3 / ha_constant
                conduction_band_unit = 'Ha'
            if conduction_band_unit != 'Ha':
                raise ValueError("Invalid conduction band unit: %s." % conduction_band_unit)

            d_2 = (y_3 - y_2) / (x_3 - x_2)
            d_1 = (y_2 - y_1) / (x_2 - x_1)
            dd = (d_2 - d_1) / ((x_3 - x_1) / 2)
            mass = 1 / dd
            elec_effective_masses.append((dir_name, mass))

        if not elec_effective_masses:
            raise ValueError("Calculating electron effective mass failed, no electron effective mass is calculated.")

        self.output('')
        self.output('(2) Electrons on conduction band.')
        self.output('Min value on conduction band: %s %s' % (elec_value, band_gap_unit))
        self.output("CBM position: %s" % elec_position)
        dir_content = ""
        for elec_dir in elec_directions:
            dir_content += "  %s(%s) -> %s(%s)\n" % (elec_dir[0][0], elec_dir[0][1], elec_dir[1][0], elec_dir[1][1])
        self.output("Electron direction(s):\n%s" % dir_content[:-1])
        mass_content = ""
        for elec_mass in elec_effective_masses:
            mass_content += "  Direction: %s, Mass: %s m0\n" % (elec_mass[0], elec_mass[1])
        self.output("Electron effective mass(es):\n%s" % mass_content[:-1])


if __name__ == '__main__':
    opts, args = getopt.getopt(sys.argv[1:], "p:s:c:hv", ["project=", "xsd=", "bs_xcd=", "help", "version"])
    project_directory = ""
    model_xsd_file = "default"
    bandstr_xcd_file = "default"

    for op, value in opts:
        if op in ("-h", "--help"):
            print(__str__)
            exit()
        elif op in ("-v", "--version"):
            print("EffectiveMass %s" % __version__)
            exit()
        elif op in ("-p", "--project"):
            project_directory = value
        elif op in ("-s", "--xsd"):
            model_xsd_file = value
        elif op in ("-c", "--xcd"):
            bandstr_xcd_file = value
        else:
            raise ValueError("Invalid option: %s." % op)

    if not project_directory:
        raise ValueError("Project directory is empty.")

    em = EffectiveMass(project_directory)
    c_paras, k_paras, b_paras = em.load(xsd_file=model_xsd_file, bs_xcd_file=bandstr_xcd_file)
    em.calculate(cell_vectors=c_paras, k_path_paras=k_paras, band_paras=b_paras)
