#include "mainwindow.h"
#include "ui_mainwindow.h"


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    curve = new QwtPlotCurve * [NUM];
    mark = new QwtPlotMarker * [NUM];

    for (int j = 0; j < NUM; ++j)
    {
        curve[j] = new QwtPlotCurve;
        curve[j]->setRenderHint(QwtPlotItem::RenderAntialiased);

        mark[j] = new QwtPlotMarker;
        mark[j]->setSymbol(new QwtSymbol(QwtSymbol::Ellipse, QColor(Qt::transparent), QColor(Qt::black), QSize(5,5)));
    }

    ui->qwtPlot->setAxisScale(QwtPlot::xBottom, 0, F_WIDTH);
    ui->qwtPlot->setAxisScale(QwtPlot::yLeft, 0, F_HEIGHT);
    grid = new QwtPlotGrid;
    grid->setMajPen(QPen(Qt::black,0,Qt::DotLine));
    grid->attach(ui->qwtPlot);

    ///////////////////////////////////////
    curve2 = new QwtPlotCurve * [NUM];
    mark2 = new QwtPlotMarker * [NUM];

    for (int j = 0; j < NUM; ++j)
    {
        curve2[j] = new QwtPlotCurve;
        curve2[j]->setRenderHint(QwtPlotItem::RenderAntialiased);

        mark2[j] = new QwtPlotMarker;
        mark2[j]->setSymbol(new QwtSymbol(QwtSymbol::Ellipse, QColor(Qt::transparent), QColor(Qt::black), QSize(5,5)));
    }

    ui->qwtPlot2->setAxisScale(QwtPlot::xBottom, -20, 20);
    ui->qwtPlot2->setAxisScale(QwtPlot::yLeft, -20, 20);
    grid2 = new QwtPlotGrid;
    grid2->setMajPen(QPen(Qt::black,0,Qt::DotLine));
    grid2->attach(ui->qwtPlot2);
    ///////////////////////////////////////

    v = 2;
    a = 0.16;
    r = 5;
    Tmin = 0;
    Tmax = 200;
    count = 100000;
    t = 0;
    connect(this, SIGNAL(stepChanged(int)), ui->lcdNumber, SLOT(display(int)));
    ui->lcdNumber->setNumDigits(trunc(log10(count)));

    w = new double[NUM];
    for (int i = 0; i < NUM; ++i)
        w[i] = (double)rand()/RAND_MAX * 0.1 + 0.95;

    x = new double * [NUM];
    y = new double * [NUM];

    xr = new double * [NUM];
    yr = new double * [NUM];
    zr = new double * [NUM];

    f = new double [NUM];

    for (int j = 0; j < NUM; ++j)
    {
        x[j] = new double[count];
        y[j] = new double[count];

        x[j][0] = (double)rand()/RAND_MAX * F_WIDTH;
        y[j][0] = (double)rand()/RAND_MAX * F_HEIGHT;

        xr[j] = new double[count];
        yr[j] = new double[count];
        zr[j] = new double[count];

        xr[j][0] = -8 + (double)rand()/RAND_MAX * 16;
        yr[j][0] = -8 + (double)rand()/RAND_MAX * 16;
        zr[j][0] = -2 + (double)rand()/RAND_MAX * 4;

        f[j] = atan(yr[j][0]/xr[j][0]);
    }

    used.resize(NUM);

    dt = (Tmax - Tmin) / 10000;  // /count;

    timer = new QTimer(this);
    connect(timer, SIGNAL(timeout()), this, SLOT(replot()));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::solveStep(int i)
{
    double X1[NUM], X2[NUM], X3[NUM], X4[NUM];
    double Y1[NUM], Y2[NUM], Y3[NUM], Y4[NUM];
    double Z1[NUM], Z2[NUM], Z3[NUM], Z4[NUM];

    double X1_[NUM], X2_[NUM], X3_[NUM];
    double Y1_[NUM], Y2_[NUM], Y3_[NUM];
    double Z1_[NUM], Z2_[NUM], Z3_[NUM];

    double x_[NUM], y_[NUM], z_[NUM];

    used.clear();
    used.resize(NUM);

    for (int j = 0; j < NUM; ++j)
    {
        x_[j] = xr[j][i-1];
        y_[j] = yr[j][i-1];
        z_[j] = zr[j][i-1];

        d[j].clear();

        for (int k = 0; k < NUM; ++k)
        {
            if ((x[k][i-1] - x[j][i-1]) * (x[k][i-1] - x[j][i-1])
                    + (y[k][i-1] - y[j][i-1]) * (y[k][i-1] - y[j][i-1]) < r * r)
                d[j].append(k);
        }
    }

    for (int j = 0; j < NUM; ++j)
    {
        X1[j] = dx(x_, y_, z_, j, i) * dt;
        Y1[j] = dy(x_, y_, z_, j, i) * dt;
        Z1[j] = dz(x_, y_, z_, j, i) * dt;

        X1_[j] = xr[j][i-1] + X1[j] / 2;
        Y1_[j] = yr[j][i-1] + Y1[j] / 2;
        Z1_[j] = zr[j][i-1] + Z1[j] / 2;
    }

    for (int j = 0; j < NUM; ++j)
    {
        X2[j] =dx(X1_, Y1_, Z1_, j, i) * dt;
        Y2[j] =dy(X1_, Y1_, Z1_, j, i) * dt;
        Z2[j] =dz(X1_, Y1_, Z1_, j, i) * dt;

        X2_[j] = xr[j][i-1] + X2[j] / 2;
        Y2_[j] = yr[j][i-1] + Y2[j] / 2;
        Z2_[j] = zr[j][i-1] + Z2[j] / 2;
    }

    for (int j = 0; j < NUM; ++j)
    {
        X3[j] = dx(X2_, Y2_, Z2_, j, i) * dt;
        Y3[j] = dy(X2_, Y2_, Z2_, j, i) * dt;
        Z3[j] = dz(X2_, Y2_, Z2_, j, i) * dt;

        X3_[j] = xr[j][i-1] + X3[j];
        Y3_[j] = yr[j][i-1] + Y3[j];
        Z3_[j] = zr[j][i-1] + Z3[j];
    }

    for (int j = 0; j < NUM; ++j)
    {
        X4[j] = dx(X3_, Y3_, Z3_, j, i) * dt;
        Y4[j] = dy(X3_, Y3_, Z3_, j, i) * dt;
        Z4[j] = dz(X3_, Y3_, Z3_, j, i) * dt;

        xr[j][i] = xr[j][i-1] + (X1[j] + 2 * X2[j] + 2 * X3[j] + X4[j]) / 6;
        yr[j][i] = yr[j][i-1] + (Y1[j] + 2 * Y2[j] + 2 * Y3[j] + Y4[j]) / 6;
        zr[j][i] = zr[j][i-1] + (Z1[j] + 2 * Z2[j] + 2 * Z3[j] + Z4[j]) / 6;

        f[j] = atan(yr[j][i]/xr[j][i]);

        if (!used.testBit(j))
        {
            comp.clear();
            dfs(j);
            double theta1 = (double)rand()/RAND_MAX * 2 * PI - PI;
//            qDebug()<<"===";
            for (int k = 0; k < comp.size(); ++k)
            {
//                qDebug()<<comp.at(k)<<f[comp.at(k)];
                for (int l = k+1; l < comp.size(); ++l)
                {
                    if (fabs(f[comp.at(k)] - f[comp.at(l)]) < 0.01)
                        qDebug()<<comp.at(k)<<"="<<comp.at(l);
                }
            }
        }

//        qDebug()<<"["<<j<<"] "<<f[j];
    }

    double *theta = new double [NUM];
    for (int j = 0; j < NUM; ++j)
    {
        theta[j] = (double)rand()/RAND_MAX * 2 * PI - PI;
        double vx = v * cos(theta[j]);
        double vy = v * sin(theta[j]);

        x[j][i] = x[j][i-1] + vx;
        y[j][i] = y[j][i-1] + vy;

        if (x[j][i] < 0)
            x[j][i] = -x[j][i];

        if (x[j][i] > F_WIDTH)
            x[j][i] = 2 * F_WIDTH - x[j][i];

        if (y[j][i] < 0)
            y[j][i] = -y[j][i];

        if (y[j][i] > F_HEIGHT)
            y[j][i] = 2 * F_HEIGHT - y[j][i];
    }
}

void MainWindow::replot()
{
    emit stepChanged(t);
    t++;
    solveStep(t);

    for (int j = 0; j < NUM; ++j)
    {
        QVector <double> xtrail, ytrail;
        int ind = std::max(0, t - 3000 / (int)Tmax);
        for (int l = ind; l < t; l++)
        {
            xtrail.append(x[j][l]);
            ytrail.append(y[j][l]);
        }

        curve[j]->setSamples(xtrail, ytrail);
        curve[j]->attach(ui->qwtPlot);
        mark[j]->setValue(x[j][t-1], y[j][t-1]);
        mark[j]->attach(ui->qwtPlot);
    }
    ui->qwtPlot->replot();

    //////////////////////////////////
    for (int j = 0; j < NUM; ++j)
    {
        QVector <double> xtrail, ytrail;
        int ind = std::max(0, t - 50000 / (int)Tmax);
        for (int l = ind; l < t; l++)
        {
            xtrail.append(xr[j][l]);
            ytrail.append(yr[j][l]);
        }

        curve2[j]->setSamples(xtrail, ytrail);
        curve2[j]->attach(ui->qwtPlot2);
        mark2[j]->setValue(xr[j][t-1], yr[j][t-1]);
        mark2[j]->attach(ui->qwtPlot2);
    }
    ui->qwtPlot2->replot();

    //////////////////////////////////

    if (t > count - 2)
    {
        killTimer(timer->timerId());
        ui->pauseBtn->setEnabled(false);
        return;
    }
}

void MainWindow::on_pauseBtn_clicked()
{
    if (ui->pauseBtn->text() == "Pause")
    {
        timer->stop();
        ui->pauseBtn->setText("Play");
    }
    else
    {
        timer->start(10);
        ui->pauseBtn->setText("Pause");
    }
}

void MainWindow::on_syncBox_clicked()
{
    if (ui->syncBox->isChecked());//
//        d = 0.5;
    else;//
//        d = 0;
}

void MainWindow::dfs(int k)
{
    used.setBit(k,true);
    comp.append(k);
    for (int i = 0; i < d[k].size(); ++i)
    {
        int to = d[k].at(i);
        if (!used.testBit(to))
            dfs(to);
    }
}
