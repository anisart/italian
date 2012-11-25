#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QKeyEvent>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    curve = new QwtPlotCurve * [num];
    mark = new QwtPlotMarker * [num];
    mark1 = new QwtPlotMarker * [num];

    int rColor[num];
    int gColor[num];
    int bColor[num];

    for (int j = 0; j < num; ++j)
    {
        rColor[j] = 50 + (double)rand()/RAND_MAX * 205;
        gColor[j] = 50 + (double)rand()/RAND_MAX * 205;
        bColor[j] = 50 + (double)rand()/RAND_MAX * 205;

        curve[j] = new QwtPlotCurve;
        curve[j]->setRenderHint(QwtPlotItem::RenderAntialiased);
        curve[j]->setPen(QColor(rColor[j],gColor[j],bColor[j]));

        mark[j] = new QwtPlotMarker;
        mark[j]->setSymbol(new QwtSymbol(QwtSymbol::Ellipse, QColor(rColor[j],gColor[j],bColor[j]), QColor(Qt::black), QSize(5,5)));
        mark1[j] = new QwtPlotMarker;
//        mark1[j]->setSymbol(new QwtSymbol());
    }

    ui->qwtPlot->setAxisScale(QwtPlot::xBottom, 0, areaWidth);
    ui->qwtPlot->setAxisScale(QwtPlot::yLeft, 0, areaHeight);
    grid = new QwtPlotGrid;
    grid->setMajPen(QPen(Qt::black,0,Qt::DotLine));
    grid->attach(ui->qwtPlot);

    ///////////////////////////////////////
    curve2 = new QwtPlotCurve * [num];
    mark2 = new QwtPlotMarker * [num];

    for (int j = 0; j < num; ++j)
    {
        curve2[j] = new QwtPlotCurve;
        curve2[j]->setRenderHint(QwtPlotItem::RenderAntialiased);
        curve2[j]->setPen(QColor(rColor[j],gColor[j],bColor[j]));

        mark2[j] = new QwtPlotMarker;
        mark2[j]->setSymbol(new QwtSymbol(QwtSymbol::Ellipse, QColor(rColor[j],gColor[j],bColor[j]), QColor(Qt::black), QSize(5,5)));
    }

    ui->qwtPlot2->setAxisScale(QwtPlot::xBottom, -30, 30);
    ui->qwtPlot2->setAxisScale(QwtPlot::yLeft, -30, 30);
    grid2 = new QwtPlotGrid;
    grid2->setMajPen(QPen(Qt::black,0,Qt::DotLine));
    grid2->attach(ui->qwtPlot2);
    ///////////////////////////////////////

    speed = ui->speedSld->value();
    connect(this, SIGNAL(stepChanged(int)), ui->lcdNumber, SLOT(display(int)));
    ui->lcdNumber->setNumDigits(trunc(log10(count)));
    connect(ui->speedSld, SIGNAL(valueChanged(int)), this, SLOT(speedChange(int)));

    used.resize(num);

    dt = (Tmax - Tmin) / (double)10000;  // /count;

    timer = new QTimer(this);
    connect(timer, SIGNAL(timeout()), this, SLOT(replot()));

    //////////////////////

    w = new double[num];

    x = new double * [num];
    y = new double * [num];

    xr = new double * [num];
    yr = new double * [num];
    zr = new double * [num];

//    f = new double [num];

    for (int j = 0; j < num; ++j)
    {
        x[j] = new double[count];
        y[j] = new double[count];

        xr[j] = new double[count];
        yr[j] = new double[count];
        zr[j] = new double[count];
    }

    flag = new TriMatrix<int>(num);
    wss = new TriMatrix<double>(num);
    on_syncBox_clicked();
    initVars();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::initVars()
{
    t = 0;

    for (int j = 0; j < num; ++j)
    {
        w[j] = (double)rand()/RAND_MAX * 0.1 + 0.95;

        x[j][0] = (double)rand()/RAND_MAX * areaWidth;
        y[j][0] = (double)rand()/RAND_MAX * areaHeight;

        xr[j][0] = -10 + (double)rand()/RAND_MAX * 20;
        yr[j][0] = -10 + (double)rand()/RAND_MAX * 20;
        zr[j][0] = 0;//-4 + (double)rand()/RAND_MAX * 8;

        ws[j] = 0;
    }

    flag->fill(0);
    wss->fill(0);

    ui->countLbl->setText(QString::number(ui->countLbl->text().toInt()+1));
}

void MainWindow::solveStep(int i)
{
    double X1[num], X2[num], X3[num], X4[num];
    double Y1[num], Y2[num], Y3[num], Y4[num];
    double Z1[num], Z2[num], Z3[num], Z4[num];

    double X1_[num], X2_[num], X3_[num];
    double Y1_[num], Y2_[num], Y3_[num];
    double Z1_[num], Z2_[num], Z3_[num];

    double x_[num], y_[num], z_[num];

    used.clear();
    used.resize(num);
    used2.clear();
    used2.resize(num);

    for (int j = 0; j < num; ++j)
    {
        x_[j] = xr[j][i-1];
        y_[j] = yr[j][i-1];
        z_[j] = zr[j][i-1];

        d[j].clear();

        for (int k = 0; k < num; ++k)
        {
            if ((x[k][i-1] - x[j][i-1]) * (x[k][i-1] - x[j][i-1])
                    + (y[k][i-1] - y[j][i-1]) * (y[k][i-1] - y[j][i-1]) < r * r)
                d[j].append(k);
        }
    }

    for (int j = 0; j < num; ++j)
    {
        X1[j] = dx(x_, y_, z_, j, i) * dt;
        Y1[j] = dy(x_, y_, z_, j, i) * dt;
        Z1[j] = dz(x_, y_, z_, j, i) * dt;

        X1_[j] = xr[j][i-1] + X1[j] / 2;
        Y1_[j] = yr[j][i-1] + Y1[j] / 2;
        Z1_[j] = zr[j][i-1] + Z1[j] / 2;
    }

    for (int j = 0; j < num; ++j)
    {
        X2[j] = dx(X1_, Y1_, Z1_, j, i) * dt;
        Y2[j] = dy(X1_, Y1_, Z1_, j, i) * dt;
        Z2[j] = dz(X1_, Y1_, Z1_, j, i) * dt;

        X2_[j] = xr[j][i-1] + X2[j] / 2;
        Y2_[j] = yr[j][i-1] + Y2[j] / 2;
        Z2_[j] = zr[j][i-1] + Z2[j] / 2;
    }

    for (int j = 0; j < num; ++j)
    {
        X3[j] = dx(X2_, Y2_, Z2_, j, i) * dt;
        Y3[j] = dy(X2_, Y2_, Z2_, j, i) * dt;
        Z3[j] = dz(X2_, Y2_, Z2_, j, i) * dt;

        X3_[j] = xr[j][i-1] + X3[j];
        Y3_[j] = yr[j][i-1] + Y3[j];
        Z3_[j] = zr[j][i-1] + Z3[j];
    }

    for (int j = 0; j < num; ++j)
    {
        X4[j] = dx(X3_, Y3_, Z3_, j, i) * dt;
        Y4[j] = dy(X3_, Y3_, Z3_, j, i) * dt;
        Z4[j] = dz(X3_, Y3_, Z3_, j, i) * dt;

        xr[j][i] = xr[j][i-1] + (X1[j] + 2 * X2[j] + 2 * X3[j] + X4[j]) / 6;
        yr[j][i] = yr[j][i-1] + (Y1[j] + 2 * Y2[j] + 2 * Y3[j] + Y4[j]) / 6;
        zr[j][i] = zr[j][i-1] + (Z1[j] + 2 * Z2[j] + 2 * Z3[j] + Z4[j]) / 6;

        double X = xr[j][i]; double Y = yr[j][i]; double Z = zr[j][i]; double W = w[j];
        ws[j] = ((Z+W*Y)*(Y*W*W-X*a*W-Y*a*a+Z*W)+(a*Y+W*X)*(X*W*W+a*Y*W+b-c*Z+X*Z))/((Z+W*Y)*(Z+W*Y)+(a*Y+W*X)*(a*Y+W*X));
//        if(i%10==0)
//        {
//            ws[j] /= 10;
////            qDebug()<<"["<<j<<"] "<<f[j]<<ws[j];
//            ws[j] = 0;
//        }
    }

    double *theta = new double [num];
    for (int j = 0; j < num; ++j)
    {
        if (!used.testBit(j))
        {
            comp.clear();
            dfs(j);
//            double theta1 = (double)rand()/RAND_MAX * 2 * PI - PI;
//            qDebug()<<"===";
            for (int k = 0; k < comp.size(); ++k)
            {
                int K = comp.at(k);
//                theta[K] = (double)rand()/RAND_MAX * 2 * PI - PI;
//                qDebug()<<K;//<<f[K];
                d2[K].clear();
                for (int l = 0; l < comp.size(); ++l)
                {
                    int L = comp.at(l);
                    if (flag->at(K,L) < tt)
                        {
                            flag->setValue(K, L, flag->at(K,L) + 1);
                            wss->setValue(K, L, wss->at(K,L) + fabs(ws[L] - ws[K]));
                        }
                    else
                        if (flag->at(K,L) > tt)
                            d2[K].append(L);
                        else
                        {
                            wss->setValue(K, L, wss->at(K,L) + fabs(ws[L] - ws[K]));
                            if (wss->at(K,L) < 0.01)
                            {
                                flag->setValue(K, L, flag->at(K,L) + 1);
                            }
                            else
                            {
                                flag->setValue(K,L,0);
                                wss->setValue(K,L,0);
                            }
                        }
                }
            }
        }

        for (int j = 0; j < num; ++j)
        {
            if (!used2.testBit(j))
            {
                comp2.clear();
                dfs2(j);
                double angle = (double)rand()/RAND_MAX * 2 * PI - PI;
                for (int k = 0; k < comp2.size(); ++k)
                {
                    theta[comp2.at(k)] = angle;
                    if (comp2.size() > 1)
                        mark1[comp2.at(k)]->setSymbol
                            (new QwtSymbol(QwtSymbol::Ellipse, QColor(Qt::transparent), QPen(QColor(Qt::red),3), QSize(8,8)));
                    else
                        mark1[comp2.at(k)]->setSymbol(new QwtSymbol());
                }
            }
        }

//        theta[j] = (double)rand()/RAND_MAX * 2 * PI - PI;
        double vx = v * cos(theta[j]);
        double vy = v * sin(theta[j]);

        x[j][i] = x[j][i-1] + vx;
        y[j][i] = y[j][i-1] + vy;

        if (x[j][i] < 0)
            x[j][i] = -x[j][i];

        if (x[j][i] > areaWidth)
            x[j][i] = 2 * areaWidth - x[j][i];

        if (y[j][i] < 0)
            y[j][i] = -y[j][i];

        if (y[j][i] > areaHeight)
            y[j][i] = 2 * areaHeight - y[j][i];
    }
}

void MainWindow::replot()
{
    emit stepChanged(t);
    t++;
    solveStep(t);
    if (t%tv!=0) return;

    for (int j = 0; j < num; ++j)
    {
        QVector <double> xtrail, ytrail;
        int ind = std::max(0, t - 6000 / Tmax);
        for (int l = ind; l < t; l++)
        {
            xtrail.append(x[j][l]);
            ytrail.append(y[j][l]);
        }

        curve[j]->setSamples(xtrail, ytrail);
        curve[j]->attach(ui->qwtPlot);
        mark[j]->setValue(x[j][t-1], y[j][t-1]);
        mark[j]->attach(ui->qwtPlot);
        mark1[j]->setValue(x[j][t-1], y[j][t-1]);
        mark1[j]->attach(ui->qwtPlot);
    }
    ui->qwtPlot->replot();

    //////////////////////////////////
    for (int j = 0; j < num; ++j)
    {
        QVector <double> xtrail, ytrail;
        int ind = std::max(0, t - 50000 / Tmax);
        for (int l = ind; l < t; l++)
        {
            if((xr[j][l] < -30)||(xr[j][l] > 30)||(yr[j][l] < -30)||(yr[j][l] > 30))
            {
//                stopProcess();
                on_pauseBtn_clicked();
//                initVars();
                qDebug()<<j<<xr[j][l]<<yr[j][l];
            }

            xtrail.append(xr[j][l]);
            ytrail.append(yr[j][l]);//qDebug()<<xtrail;
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
        stopProcess();
//        return;
    }

//    if (t==500) initVars();
}

void MainWindow::speedChange(int value)
{
    speed = value;
    if (ui->pauseBtn->text() == "Pause")
    {
        timer->start(speed);
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
        timer->start(speed);
        ui->pauseBtn->setText("Pause");
    }
}

void MainWindow::on_syncBox_clicked()
{
    if (ui->syncBox->isChecked())
        dd = D_VAL;
    else
        dd = 0;
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

void MainWindow::dfs2(int k)
{
    used2.setBit(k,true);
    comp2.append(k);
    for (int i = 0; i < d2[k].size(); ++i)
    {
        int to = d2[k].at(i);
        if (!used2.testBit(to))
            dfs2(to);
    }
}

void MainWindow::stopProcess()
{
    killTimer(timer->timerId());
    ui->pauseBtn->setEnabled(false);
    ui->speedSld->setEnabled(false);
}

void MainWindow::keyPressEvent(QKeyEvent * event)
{
    if (event->key() == Qt::Key_Escape)
        qApp->quit();
}

void MainWindow::on_resetBtn_clicked()
{
    initVars();
    ui->pauseBtn->setText("Pause");
    ui->pauseBtn->setEnabled(true);
    ui->speedSld->setEnabled(true);
    timer->start(speed);
}
