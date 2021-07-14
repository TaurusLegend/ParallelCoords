"""
Version: 
    1.0

Copyright Notice:
    Copyright (c) 2021, Claudia Campos
    All rights reserved.

License Notice:

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""

from matplotlib import pyplot as plt
from matplotlib.backend_bases import MouseButton
from matplotlib.backend_tools import Cursors
import numpy as np


class ParallelCoordinates:
    def __init__(self, data):  # alocacao no construtor
        """
        [Description]
        
        Parameters
        ----------
        data : 2D array
        
        """
        self.data = np.array(data, dtype=np.float64)  # Dados multidimensionais - ordem original
        self.dataNorm = np.zeros(self.data.shape)  # Dados normalizados - Ordem de visualizacao
        self.nPoints, self.nDims, = self.data.shape  # Numero de pontos e numero de dimensoes
        self.maxDims = np.max(self.data, 0)  # Lista com os valores maximos para cada dimensao - ordem original
        self.minDims = np.min(self.data, 0)  # Lista com os valores minimos para cada dimensao - ordem original
        self.px = np.arange(self.nDims)  # Permutacao dos eixos
        self.init = np.arange(self.nDims)
        self.posAxis = np.arange(self.nDims, dtype=np.float64)  # Posicao de cada eixo - Ordem de visualizacao
        self.normalizeData()  # Dados normalizados - Ordem de visualizacao
        self.autoSpace = True  # Espacamento automatico - True ou False
        self.scale = 1.1  # Valor da escala
        self.defaultScale = 1.1
        self.ind = None  # Indice do eixo selecionado
        self.yPrevious = None  # Coordenada y anterior
        self.mouseX = None  # Coordenada x do rato
        self.mouseY = None  # Coordenada y do rato
        self.lines = list()  # np.zeros(self.nPoints, dtype=object)  # Lista de polilinhas (Line2D)
        self.defaultAlpha = 1
        self.alpha = 1  # Transparencia de todas as linhas do grafico
        self.color = 'b'  # Cor de todas as linhas do grafico
        self.defaultColor = 'b'
        self.fig, self.host = plt.subplots()  # Figura atual do grafico, Axes principal
        self.axis = np.zeros(self.nDims, dtype=object)  # Lista de eixos (Axes)
        self.setGraph()  # Faz setup do grafico e desenha-o

    def getColor(self):
        """
        Returns the current color of all polylines in the plot.
        
        Returns
        -------
        color : color
        
        """
        return self.color

    def setColor(self, color):
        """
        Sets the color of all polylines in the plot.
        
        Parameters
        ----------
        color : color
        
        """
        self.color = color

    def updateColor(self, color):
        """
        Updates the color of all polylines in the plot.
        
        Parameters
        ----------
        color : color
        
        """        
        self.setColor(color)
        for line in self.lines:
            line._color = self.color
        self.host.stale = True

    def getAlpha(self):
        """
        Returns the current alpha value of all polylines in the plot.
        
        Returns
        -------
        alpha : float
        
        """
        return self.alpha

    def setAlpha(self, alpha):
        """
        Sets the alpha value of all polylines in the plot.
        
        Parameters
        ----------
        alpha : float
        
        """        
        self.alpha = alpha

    def updateAlpha(self, alpha):
        """
        Updates the alpha value of all polylines in the plot.
        
        Parameters
        ----------
        alpha : float
        
        """
        self.setAlpha(alpha)
        for line in self.lines:
            line._alpha = self.alpha
        self.host.stale = True

    def getAutoSpace(self):
        """
        Returns whether the axes are automatically spaced or not.
        
        Returns
        -------
        autoSpace : bool
        
        """        
        return self.autoSpace

    def setAutoSpace(self, autoSpace):
        """
        Sets whether the axes are automatically spaced or not.
        
        Parameters
        ----------
        autoSpace : bool
        
        """           
        self.autoSpace = autoSpace

    def normalizeData(self, i=slice(None)):  # i = slice(0, self.nDims)
        """
        Normalizes the data according to the minimum and maximum for each objective.
        
        Parameters
        ----------
        i : int or slice, default: slice(None)
            Optional
        
        """
        self.dataNorm[:, i] = (self.data[:, self.px[i]] - self.minDims[self.px[i]]) / (
                self.maxDims[self.px[i]] - self.minDims[self.px[i]])

    def updateLinesDataXY(self):
        """
        Invalidates the x and y data for each line present in the plot and forces the plot to be redraw again.
        """
        for line in self.lines:
            line._invalidy = True
            line._invalidx = True

        self.host.stale = True
        plt.draw()

    def updateLinesDataY(self):
        """
        Invalidates the y data for each line present in the plot and forces the plot to be redraw again.
        """
        for line in self.lines:
            line._invalidy = True

        self.host.stale = True
        plt.draw()

    def removeLines(self):
        """
        Remove all the lines from the plot without removing the axes or closing the figure.
        """
        for line in self.lines:
            line.remove()

        self.lines.clear()

    # Adaptado de: https://matplotlib.org/2.0.2/examples/pylab_examples/multiple_yaxis_with_spines.html
    def makePatchSpinesInvisible(self, axes):
        """
        Make patch and all spines invisible for the giving axes.
        
        Parameters
        ----------
        axes : Axes
        
        """
        axes.set_frame_on(True)
        axes.patch.set_visible(False)
        for spine in axes.spines.values():
            spine.set_visible(False)

    def setData(self, data):
        """
        Replace the data in the current figure without changing the number of objectives.
        
        Parameters
        ----------
        data : 2D array
        
        """
        if self.data.shape == data.shape:
            self.data[:] = data[:]
            self.maxDims = np.max(self.data, 0)
            self.minDims = np.min(self.data, 0)
            self.setPermutation(self.init)
            self.setPosition(self.init)
            self.normalizeData()

            for px, pos in zip(self.px, self.posAxis):
                self.axis[px].spines["right"].set_position(("data", pos))
                self.axis[px].set_ylim(self.minDims[px], self.maxDims[px])

            self.updateLinesDataXY()

        elif data.shape[1] == self.nDims:
            self.data = data
            self.dataNorm = np.zeros(self.data.shape)
            self.nPoints, _, = self.data.shape
            self.maxDims = np.max(self.data, 0)
            self.minDims = np.min(self.data, 0)
            self.setPermutation(self.init)
            self.setPosition(self.init)
            self.normalizeData()

            for px, pos in zip(self.px, self.posAxis):
                self.axis[px].spines["right"].set_position(("data", pos))
                self.axis[px].set_ylim(self.minDims[px], self.maxDims[px])

            self.removeLines()

            self.lines = self.host.plot(self.posAxis, self.dataNorm.transpose(), color=self.color, alpha=self.alpha)

            plt.draw()

        else:
            raise RuntimeError()

    def setPx(self, px, i=slice(None)):
        """
        Sets the permutation vector that swaps the parallel axes and updates the plot.
        
        Parameters
        ----------
        px : 1D array
        
        i : int or slice, default: slice(None)
            Optional
        
        """
        self.px[i] = px

    def setPosAxis(self, posAxis, i=slice(None)):
        """
        Sets the position vector.
        
        Parameters
        ----------
        px : 1D array
        
        i : int or slice, default: slice(None)
            Optional
        
        """
        self.posAxis[i] = posAxis

    def setPermutation(self, px):
        """
        Sets the permutation array that swaps the parallel axes and updates the plot.
        
        Parameters
        ----------
        px : 1D array
        
        Notes
        -----
        (...)
        
        """        
        self.setPx(px)

        for px, pos in zip(self.px, self.posAxis):
            self.axis[px].spines["right"].set_position(("data", pos))

        self.normalizeData()
        self.updateLinesDataXY()

    def setPosition(self, posAxis):
        """
        Sets the position for each axis in visualization order.
        
        Parameters
        ----------
        posAxis : 1D array
        
        Notes
        -----
        TODO Explain better this concept
        
        """
        self.setPx(self.px[np.argsort(posAxis)])

        if self.autoSpace:
            self.setPosAxis(self.init)
        else:
            self.setPosAxis(np.sort(posAxis))

        for px, pos in zip(self.px, self.posAxis):
            self.axis[px].spines["right"].set_position(("data", pos))

        self.normalizeData()
        self.updateLinesDataXY()

    def getScale(self):
        """
        Returns the current zoom scale.
        
        Returns
        -------
        scale : float
        
        """
        return self.scale

    def setScale(self, scale):
        """
        Sets the zoom scale.
        
        Parameters
        ----------
        scale : float
        
        """
        self.scale = scale

    def resetAttributes(self):
        """
        Resets all attributes to all their default values.
        """
        self.setPx(self.init)
        self.setPosAxis(self.init)
        
        self.maxDims[:] = np.max(self.data, 0)
        self.minDims[:] = np.min(self.data, 0)      
        
        if self.color != self.defaultColor:
            self.setColor(self.defaultColor)
        if self.alpha != self.defaultAlpha:
            self.setAlpha(self.defaultAlpha)
        if self.scale != self.defaultScale:
            self.setScale(self.defaultScale)
        if not self.autoSpace:
            self.setAutoSpace(True)

        for px, pos in zip(self.px, self.posAxis):
            self.axis[px].spines["right"].set_position(("data", pos))
            self.axis[px].set_ylim(self.minDims[px], self.maxDims[px])

        self.normalizeData()
        
        for line in self.lines:
            line._invalidy = True
            line._invalidx = True
            line._alpha = self.alpha
            line._color = self.color
            
        self.host.stale = True
        plt.draw()

    def setGraph(self):
        """
        Sets the axes in the correct position with the correct attributes, makes the plot and attach the listeners to the events.
        """
        self.makePatchSpinesInvisible(self.host)
        self.host.get_xaxis().set_ticks(list())
        self.host.axes.get_xaxis().set_visible(False)
        self.host.axes.get_yaxis().set_visible(False)
        self.host.set_xlim(- 0.5, self.nDims - 0.5)
        self.host.set_ylim(0, 1)

        for i in range(self.nDims):
            self.axis[i] = self.host.twinx()
            self.makePatchSpinesInvisible(self.axis[i])
            self.axis[i].spines["right"].set_visible(True)
            self.axis[i].spines["right"].set_position(("data", self.posAxis[i]))
            self.axis[i].set_ylim(self.minDims[i], self.maxDims[i])

        self.lines = self.host.plot(self.posAxis, self.dataNorm.transpose(), color=self.color, alpha=self.alpha)

        plt.ion()
        self.fig.canvas.mpl_connect('button_press_event', self.pressMouseButton)
        self.fig.canvas.mpl_connect('button_release_event', self.releaseMouseButton)
        self.fig.canvas.mpl_connect('motion_notify_event', self.onMouseMove)
        self.fig.canvas.mpl_connect('scroll_event', self.onScroll)
        plt.show()

    def getIndUnderAxis(self, epsilon=0.1):
        """
        Get the index of the axes that is currently being manipulated.
        
        Parameters
        ----------
        epsilon : float, default: 0.1
        
        """
        diff = np.abs(self.posAxis - self.mouseX)
        diff = np.logical_and(diff == min(diff), diff < epsilon)
        result, = np.where(diff)  # np.where(diff == True)
        if result.size == 0:
            self.ind = None
        else:
            self.ind = result[0]

    def pressMouseButton(self, event):
        """
        
        """
        if event.inaxes is None or event.button != MouseButton.LEFT:
            return
        self.mouseX, self.mouseY = self.host.transData.inverted().transform([event.x, event.y])
        self.yPrevious = self.mouseY
        self.getIndUnderAxis()

    def releaseMouseButton(self, event):
        """
        
        """
        if event.button == MouseButton.LEFT:
            if self.ind is not None:
                pos = self.posAxis[self.ind]
                if self.ind > 0 and self.posAxis[self.ind - 1] > self.posAxis[self.ind]:
                    indexes, = np.where(self.posAxis < self.posAxis[self.ind])
                    if indexes.size == 0:
                        newInd = 0
                    else:
                        newInd = indexes[indexes.size - 1] + 1
                    temp = self.px[self.ind]
                    self.px[newInd + 1:self.ind + 1] = self.px[newInd:self.ind]
                    self.px[newInd] = temp
                    if self.autoSpace:
                        self.posAxis[newInd + 1:self.ind + 1] = self.posAxis[newInd:self.ind] + 1
                        self.posAxis[newInd] = newInd
                        for i in range(newInd, self.ind + 1):
                            self.axis[self.px[i]].spines["right"].set_position(("data", self.posAxis[i]))
                    else:
                        self.posAxis[newInd + 1:self.ind + 1] = self.posAxis[newInd:self.ind]
                        self.posAxis[newInd] = pos
                    self.ind = newInd

                elif self.ind < self.nDims - 1 and self.posAxis[self.ind] > self.posAxis[self.ind + 1]:
                    indexes, = np.where(self.posAxis > self.posAxis[self.ind])
                    if indexes.size == 0:
                        newInd = self.nDims - 1
                    else:
                        newInd = indexes[0] - 1
                    temp = self.px[self.ind]
                    self.px[self.ind:newInd] = self.px[self.ind + 1:newInd + 1]
                    self.px[newInd] = temp
                    if self.autoSpace:
                        self.posAxis[self.ind:newInd] = self.posAxis[self.ind + 1:newInd + 1] - 1
                        self.posAxis[newInd] = newInd
                        for i in range(self.ind - 1, newInd + 1):
                            self.axis[self.px[i]].spines["right"].set_position(("data", self.posAxis[i]))
                    else:
                        self.posAxis[self.ind:newInd] = self.posAxis[self.ind + 1:newInd + 1]
                        self.posAxis[newInd] = pos
                    self.ind = newInd
                else:
                    if self.autoSpace:
                        self.posAxis[self.ind] = self.ind
                        self.axis[self.px[self.ind]].spines["right"].set_position(("data", self.posAxis[self.ind]))

            self.normalizeData()
            self.updateLinesDataXY()
            self.ind = None

    def onMouseMove(self, event):
        """
        
        """
        if self.ind is None:
            """for axis in self.axis:
                help(axis)
                #if axis.get_yticks().contains(event):
                    #self.fig.canvas.toolbar.set_cursor(Cursors.MOVE)
            self.fig.canvas.toolbar.set_cursor(Cursors.POINTER)"""
        elif event.button == MouseButton.LEFT:
            self.mouseX, self.mouseY = self.host.transData.inverted().transform([event.x, event.y])
            self.posAxis[self.ind] = self.mouseX
            px = self.px[self.ind]
            self.axis[px].spines["right"].set_position(("data", self.posAxis[self.ind]))

            deltaY = self.mouseY - self.yPrevious

            diff = deltaY * (self.maxDims[px] - self.minDims[px])
            self.minDims[px] -= diff
            self.maxDims[px] -= diff
            self.axis[px].set_ylim(self.minDims[px], self.maxDims[px])
            self.normalizeData(self.ind)
            self.yPrevious = self.mouseY
            self.updateLinesDataXY()

    def onScroll(self, event):
        """
        
        """
        self.mouseX, self.mouseY = self.host.transData.inverted().transform([event.x, event.y])
        self.ind = int(round(self.mouseX))

        if self.ind is not None:
            if event.button == 'up':
                scaleFactor = self.scale
            elif event.button == 'down':
                scaleFactor = 1 / self.scale
            else:
                raise RuntimeError()  # Nunca deve ser alcancado

            px = self.px[self.ind]
            yy = self.minDims[px] + self.mouseY * (self.maxDims[px] - self.minDims[px])
            self.minDims[px] = yy + (self.minDims[px] - yy) / scaleFactor
            self.maxDims[px] = yy + (self.maxDims[px] - yy) / scaleFactor
            self.axis[px].set_ylim(self.minDims[px], self.maxDims[px])

            self.normalizeData(self.ind)

            self.updateLinesDataY()
