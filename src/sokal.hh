/*
 *  This file is part of BGX, the Bayesian Gene eXpression program.
 *  Copyright (c) 2003-2004  Graeme Ambler <graeme@ambler.me.uk>
 *
 *  BGX is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License, version 2, as
 *  published by the Free Software Foundation.
 *
 *  BGX is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */

// Estimate integrated autocorrelation time using the method of Sokal.
// Note that this is actually the sum from -infty to +infty of the acf, hence 
// twice Sokal's definition.  The implementation follows that of 
// Green and Han (1990).
// WARNING: The input array x is destroyed in the process.
#ifndef _SOKAL_H
#define _SOKAL_H

extern "C" {
int sokal(int *n, double *x, double *var, double *tau, int* m);
}
#endif // _SOKAL_H
