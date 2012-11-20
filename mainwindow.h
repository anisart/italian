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

#define F_WIDTH     100
#define F_HEIGHT    100
#define PI          3.14159265
#define NUM         10
#define D_VAL       2          ///////

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
    QwtPlotMarker **mark, **mark2;
    QTimer *timer;
    QwtPlotGrid *grid, *grid2;
//    double d[NUM][NUM];
    double ws[NUM];
    int counter[NUM][NUM], flag[NUM][NUM];

    double **x, **y, **xr, **yr, **zr;
    double *w, *f;
    double a, b, c, r, v;
    double Tmin, Tmax, dt;
    int count;
    int t, speed;

    QBitArray used, used2;
    QVector <int> comp, comp2;
    QVector <int> d[NUM], d2[NUM];

    void initVars();
    void solveStep(int i);
    void dfs(int k);
    void dfs2(int k);
    void stopProcess();

    double dx(double *X, double *Y, double *Z, int k, int l)
    {
            return -w[k] * Y[k] - Z[k];
    }


    double dy(double *X, double *Y, double *Z, int k, int l)
    {
            double sync = 0;
            int i;
            foreach (i, d[k])
                sync += D_VAL * (Y[i] - Y[k]);

            return w[k] * X[k] + a * Y[k] + sync;
    }

    double dz(double *X, double *Y, double *Z, int k, int l)
    {
            return b + Z[k] * (X[k] - c);
    }

    friend QDebug operator<< (QDebug dbg, const QBitArray& array);

public slots:
    void replot();
    void speedChange(int value);

signals:
    void stepChanged(int value);

private slots:
    void on_pauseBtn_clicked();
    void on_syncBox_clicked();

    void on_resetBtn_clicked();

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
