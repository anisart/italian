#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_marker.h>
#include <qwt_symbol.h>
#include <qwt_plot_grid.h>
#include <QBitArray>
#include "trimatrix.h"

#define         PI              3.14159265
const int       areaWidth =     100;
const int       areaHeight =    100;
const int       num =           100;
#define         D_VAL           1
const double    v =             0.2;
const double    a =             0.16;
const double    b =             0.1;
const double    c =             8.5;
const double    r =             10;
const int       Tmin =          0;
const int       Tmax =          200;
const int       count =         1000000;
const int       tt =            10;
const int       tv =            10;
const int       ttt =           100;
const double    eps =           0.02;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    
private:
    Ui::MainWindow *ui;
    QwtPlotCurve **curve, **curve2;
    QwtPlotMarker **mark, **mark2, **mark1;
    QTimer *timer;
    QwtPlotGrid *grid, *grid2;
    TriMatrix <int> *flag;

    /////
    QwtPlot *qwtPlot3; QwtPlotCurve **curve3; QwtPlotGrid *grid3; QList <double> xAxis;
    ////

    QList <double> *x, *y, *xr, *yr, *zr, *ws;
    double *w, *vx, *vy;
    double dd;
    double dt;
    int t, tick, *p;

    QBitArray used, used2;
    QVector <int> comp, comp2;
    QVector <int> d[num], d2[num];

    void initVars();
    void solveStep(int i);
    void dfs(int k);
    void dfs2(int k);
    void stopProcess();

    double dx(double *X, double *Y, double *Z, int k, int l)
    {
        double sync = 0;
        foreach (int i, d[k])
            sync += dd * (X[i] - X[k]);
        return -w[k] * Y[k] - Z[k] + sync;
    }


    double dy(double *X, double *Y, double *Z, int k, int l)
    {
        double sync = 0;
        foreach (int i, d[k])
            sync += dd * (Y[i] - Y[k]);

        return w[k] * X[k] + a * Y[k] + sync;
    }

    double dz(double *X, double *Y, double *Z, int k, int l)
    {
        return b + Z[k] * (X[k] - c);
    }

    friend QDebug operator<< (QDebug dbg, const QBitArray& array);

signals:
    void stepChanged(int value);

private slots:
    void replot();
    void on_pauseBtn_clicked();
    void on_resetBtn_clicked();
    void on_forceSld_valueChanged(int value);
    void on_tickSld_valueChanged(int value);

protected:
    void keyPressEvent(QKeyEvent * event);
};

inline QDebug operator<< (QDebug dbg, const QBitArray& array)
{
    QString text;
    for (int i = 0; i < array.size(); ++i)
        text += array.testBit(i) ? "1": "0";
    dbg << text;
    return dbg;
}

#endif // MAINWINDOW_H
