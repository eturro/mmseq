#!/bin/bash
##    Copyright 2010 Ernest Turro
##
##    This program is free software; you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program; if not, write to the Free Software
##    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

ff=1
while getopts "v" flag; do
  case "$flag" in
  v) ff=0;;
  ?) printf "Usage: %s: [-v] pattern fasta_file\n" $(basename $0) >&2; exit 1;;
  esac
done

shift $(($OPTIND - 1))

awk -v pattern=$1 -v ff=$ff '
BEGIN { f=0 }
/^>/  { $0 ~ pattern ? f=1 : f=0 }
f==ff{ print }' $2
