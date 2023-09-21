//
//  main.c
//  Stable Fluids 3D
//

#include "visualizer.h"


Visualizer::Visualizer()
{
    windowWidth = 600;
    windowHeight = 600;

    translationX;
    translationY;
    translationZ;

    rotationX;
    rotationY;

    maxN = nx > ny ? (nx > nz ? nx : nz) : (ny > nz ? ny : nz);

    pause = false;
    pauseflag = true;
    debug = false;


    alpha = 0.02f;
    drawMode = 0;
    draw = true;
    drawPerspective = true;
    drawAxis = true;
    drawObstacle = true;

    logFrameStep = 0;

    total = 0.f;
    oldTotal = 0.f;

    drawScale = 0.2f;
    drawThick = 20;
    drawVolStart = (float)(nz / 2 - drawThick / 2);

    preAvgTemp = 0;

    simlog = false;
    simlogname[512];

    drawY = 1;

    room = 0;

    particleSpeed;

    temAvg = 0.f;
    ventDir = 1;
}
Visualizer::~Visualizer()
{

}

void Visualizer::initParticles() {
    srand(time(NULL));

    particleSpeed = 1.f;
}

string Visualizer::getTime()
{
    auto now = std::chrono::system_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    tm* tstruct = new tm();
    localtime_s(tstruct, &in_time_t);

    std::stringstream ss;
    //ss << std::put_time(tstruct, "%Y-%m-%d %X");
    //ss << ":" << std::setw(3) << std::setfill('0') << ms.count();
    ss << std::put_time(tstruct, "%X") << std::setw(3) << std::setfill('0') << ms.count();

    delete tstruct;
    return ss.str();
}

float Visualizer::norm2tem(float t)
{
    //EulerianSetting& eu = Setting::Instance()->eulerianSetting;
    //return (1 - t) * (eu.maxTem - eu.minTem) + eu.minTem;
    return t;
}

float Visualizer::getLocalTem()
{
    double temp = 0.f;

#pragma omp parallel for reduction (+:temp)
    for (int i = -5; i <= 5; i++)
    {
        for (int j = -5; j <= 5; j++)
        {
            for (int k = -5; k <= 5; k++)
            {
                temp += solver->getDensityGPU(px + i, py + j, pz + k);
            }
        }
    }

    
    return norm2tem(temp / (11 * 11 * 11));
}

float Visualizer::getAvgTem()
{   //room==2
    if (room == 1)
    {
        double temp = 0.f;

        // for 2nd test room
        float a = -4.78571f;
        float b = 72.78571f;

        float a2 = 6.0909f;
        float b2 = -413.182f;

        //float temp = 0.f;
        int aaa = 0;

        for (int x = 1; x < nx; x++)
        {
            //int minz = (x - b) / a;
            //int maxz = (x - b2) / a2;
            for (int z = 1; z < nz; z++)
            {
                //if (minz < z && z < maxz)
                {
                    for (int y = 1; y < ny; y++)
                    {

                        temp += solver->getDensityGPU(x, y, z);
                        aaa++;
                    }
                }
            }
        }

        return norm2tem(temp / aaa);
    }
    else
    {
        double temp = 0.f;
#pragma omp parallel for reduction (+:temp)
        for (int i = 1; i <= nx - 2; i++)
        {
            for (int j = 1; j <= ny - 2; j++)
            {
                for (int k = 1; k <= nz - 2; k++)
                {
                    temp += solver->getDensityGPU(i, j, k);
                }
            }
        }

        return norm2tem(temp / ((nx - 2) * (ny - 2) * (nz - 2)));
    }
}

void Visualizer::saveData(char* fn)
{
    FILE* f;

    if (fopen_s(&f, fn, "w") != 0)
        return;

    for (int i = 1; i <= nx; i++)
    {
        for (int j = 1; j <= ny; j++)
        {
            for (int k = 1; k <= nz; k++)
            {
                auto t = solver->getDensityGPU(i, j, k);
                auto xx = solver->getVelocityGPUU(i, j, k);
                auto yy = solver->getVelocityGPUV(i, j, k);
                auto zz = solver->getVelocityGPUW(i, j, k);

                fprintf(f, "%f %f %f %f\n", xx, yy, zz, norm2tem(t));
            }
        }
    }

    fclose(f);
}

void Visualizer::writeParam(FILE* f)
{
    Setting& s = *(setting);
    /*
    fprintf(f, "grid : %zu x %zu x %zu\n", s.eulerianSetting.gridX, s.eulerianSetting.gridY, s.eulerianSetting.gridZ);
    fprintf(f, "timeStep : %f\n", s.eulerianSetting.timeStep);
    fprintf(f, "diffuse : %f\n", s.eulerianSetting.diffuseCoef);
    fprintf(f, "viscosity : %f\n", s.eulerianSetting.viscosityCoef);
    fprintf(f, "buoyancy : %f\n", s.eulerianSetting.buoyancyCoef);
    fprintf(f, "gravity : %f\n", s.eulerianSetting.gravity);
    fprintf(f, "force : %f\n", s.eulerianSetting.force);
    fprintf(f, "srcTem : %f\n", s.eulerianSetting.srcTem);
    fprintf(f, "curTem : %f\n", s.eulerianSetting.curTem);
    fprintf(f, "tempRange : %f ~ %f\n\n", s.eulerianSetting.minTem, s.eulerianSetting.maxTem);

    fprintf(f, "machine_target : %f\n", s.machineSetting.tarTem);
    fprintf(f, "machine_ventDir : %d\n", s.machineSetting.ventDir);
    fprintf(f, "machine_ventSpeed : %d\n", s.machineSetting.ventSpeed);
    fprintf(f, "machine_initAngle : %3.2f\n", s.machineSetting.initAngle);
    fprintf(f, "machine_interval : %3.2f\n\n", s.machineSetting.interval);

    fprintf(f, "sAttenuationP : %3.2f\n", s.envSetting.sAttenuationP);
    fprintf(f, "sAttenuationR : %3.2f\n", s.envSetting.sAttenuationR);
    fprintf(f, "useUpSampling : %s\n\n", (s.envSetting.useUpSampling ? "true" : "false"));
    */
}

bool Visualizer::OpenFile(char* fn)
{
    if (!(setting->envSetting.saveData))
        return false;;

    FILE* f;
    if (fopen_s(&f, fn, "w") != 0)
        return false;

    if (setting->envSetting.saveAvg)
    {
        writeParam(f);
    }
    else
    {
        fprintf(f, "시간");
        for (int i = 0; i < 18; i++)
        {
            for (int j = 0; j < 11; j++)
            {
                for (int k = 0; k < 5; k++)
                {
                    fprintf(f, " (%d,%d,%d)", i + 1, j + 1, k + 1);
                }
            }
        }
        fprintf(f, "\n");

        // save param with apart file.
        FILE* fset;
        char setfn[100] = { 0, };
        sprintf(setfn, "set_%s", fn);

        if (fopen_s(&fset, setfn, "w") == 0)
        {
            writeParam(fset);
            fclose(fset);
        }

    }
    fclose(f);

    return true;
}

bool Visualizer::AppendData(char* fn, float time, float scale)
{
    if (!(setting->envSetting.saveData))
        return false;;

    FILE* f;

    if (fopen_s(&f, fn, "a") != 0)
        return false;

    if (setting->envSetting.saveAvg)
    {
        fprintf(f, "%f\n", preAvgTemp);
        fclose(f);
        return true;
    }

    fprintf(f, "%d", (int)(time));

    float initW = 0.8f;
    float initH = 0.1f;

    float deltaW = 0.6f;
    float deltaH = 0.5f;

    int cnt = 0;
    float total = 0.f;

    for (int i = 0; i < 18; i++)
    {
        for (int j = 0; j < 11; j++)
        {
            for (int k = 0; k < 5; k++)
            {
                int x = (int)((initW + i * deltaW) * scale);
                int z = (int)((initW + j * deltaW) * scale);

                int y = (int)((initH + k * deltaH) * scale);

                auto t = norm2tem(solver->getDensityGPU(x, y, z));
                fprintf(f, " %f", t);

                total += t;
                cnt++;
            }
        }
    }
    fprintf(f, "\n");

    fclose(f);

    //auto t = total / cnt;
    //printf("--the average temperature : %3.5f  / dif : %3.5f\n", t, t - temAvg);
    //temAvg = t;

    return true;
}

void Visualizer::help()
{
    solver->getGridSize(nx, ny, nz);

    printf("\nThe simulation size is %zd x %zd x %zd\n", nx, ny, nz);
    printf("\t a,d : rotation by Y axis\n");
    printf("\t r : reload script\n");
    printf("\t w,s : rotation by X axis\n");
    printf("\t e : reset rotation\n");

    printf("\t b : draw or not switch\n");
    printf("\t x : draw obstacle(occupancy)\n");
    printf("\t l : draw boundary(axis)\n");
    printf("\t v : draw mode switch\n");
    printf("\t 8,9 : draw scale factor adjustment\n");
    //printf("\t l : air outlet swing mode switch\n");

    printf("\t u : forcely start flow simulation\n");
    printf("\t c : reset simulation\n");
    printf("\t k : accelerator simulation\n");
    printf("\t up/down : vent angle\n");
    printf("\t ri/left : vent speed\n");
    printf("\t 1,2,3 : ac type (wall, ceiling, stand\n");
    printf("\t o : single / multiple outlet switch\n");

    printf("\t p : temperature statistic\n");

    //printf("\t 7,8,9 : timestep adjustment\n");
}

void Visualizer::wuInitialize()
{

    if (!FPP) {

        rotationX = 30;
        rotationY = -45;

        translationX = 0;
        translationY = -0.0;
        translationZ = -1.0;
    }
    else {
        rotationX = -30;
        rotationY = 0;

        translationX = 0.0;
        translationY = -0.3;
        translationZ = -0.4;
    }

    //server.setSolver(solver);
    //solver->init(nx, ny,nz);
    //server.reAllocate(98, 28, 90);

    initParticles();
}

void Visualizer::wuDrawGrid()
{
    glLineWidth(1.0f);

    float size = 1.3f;
    float x = size / 2.f;

    //glBegin(GL_QUADS);
    //glColor4f(1.0f, 1.0f, 1.0f, 0.2f);
    //glVertex3f(-x, x, -x);
    //glVertex3f(-x, x, x);
    //glVertex3f(x, x, x);
    //glVertex3f(x, x, -x);
    //glEnd();

    glBegin(GL_LINES);
    glColor3f(1.0f, 0.0f, 0.0f);

    glVertex3f(-x, 0.0f, -x);
    glVertex3f(x, 0.0f, -x);

    glColor3f(0.0f, 1.0f, 0.0f);

    glVertex3f(-x, 0.0f, -x);
    glVertex3f(-x, size, -x);

    glColor3f(0.0f, 0.0f, 1.0f);

    glVertex3f(-x, 0.0f, -x);
    glVertex3f(-x, 0.0f, x);

    glColor3f(1.0f, 0.0f, 1.0f);

    glVertex3f(x, 0.0f, -x);
    glVertex3f(x, size, -x);

    glVertex3f(x, size, -x);
    glVertex3f(-x, size, -x);

    glVertex3f(-x, size, x);
    glVertex3f(-x, 0.0f, x);

    glVertex3f(-x, size, x);
    glVertex3f(-x, size, -x);

    glVertex3f(x, 0.0f, -x);
    glVertex3f(x, 0.0f, x);

    glVertex3f(-x, 0.0f, x);
    glVertex3f(x, 0.0f, x);

    glVertex3f(x, size, -x);
    glVertex3f(x, size, x);

    glVertex3f(x, size, x);
    glVertex3f(x, 0.0f, x);

    glVertex3f(-x, size, x);
    glVertex3f(x, size, x);

    glEnd();
    /*
    glBegin(GL_LINES);
    glColor3f(1.0f, 1.0f, 1.0f);

    glVertex3f(-x, x, x);
    glVertex3f(x, x, x);

    glVertex3f(x, x, x);
    glVertex3f(x, x, -x);

    glVertex3f(x, x, -x);
    glVertex3f(-x, x, -x);

    glVertex3f(-x, x, -x);
    glVertex3f(-x, x, x);

    glVertex3f(0, x, -x);
    glVertex3f(0, x, x);

    glVertex3f(-x, x, 0);
    glVertex3f(x, x, 0);

    glEnd();
    */
}

void Visualizer::DrawObstacle()
{
    GLfloat positionX;
    GLfloat positionY;
    GLfloat positionZ;

    GLfloat scale = 1.3f;
    GLfloat h = scale / maxN;

    glBegin(GL_QUADS);
    glColor4f(1.f, 0.82f, 0.f, 0.1f);

    for (int x = 1; x <= nx + 1; x++)
    {
        positionX = (x - 1.f) * h - scale / 2.f;

        for (int y = 1; y <= ny + 1; y++)
        {
            positionY = (y - 1.f) * h;

            for (int z = 1; z <= nz + 1; z++)
            {
                positionZ = (z - 1.f) * h - scale / 2.f;

                //if (solver->getObstacleU(x, y, z)) {
                //    glVertex3f(positionX, positionY, positionZ);
                //    glVertex3f(positionX, positionY + h, positionZ);
                //    glVertex3f(positionX, positionY + h, positionZ + h);
                //    glVertex3f(positionX, positionY, positionZ + h);
                //}
                //if (solver->getObstacleV(x, y, z)) {
                //    glVertex3f(positionX, positionY, positionZ);
                //    glVertex3f(positionX + h, positionY, positionZ);
                //    glVertex3f(positionX + h, positionY, positionZ + h);
                //    glVertex3f(positionX, positionY, positionZ + h);
                //}
                //if (solver->getObstacleW(x, y, z)) {
                //    glVertex3f(positionX, positionY, positionZ);
                //    glVertex3f(positionX, positionY + h, positionZ);
                //    glVertex3f(positionX + h, positionY + h, positionZ);
                //    glVertex3f(positionX + h, positionY, positionZ);
                //}

                if (solver->getObstacle(x, y, z))
                {
                    if (drawObstacle)
                    {
                        //density000 = density010 = density100 = density110 =
                        //    density001 = density011 = density101 = density111 = 1.f;

                        //glColor4f(density111, 0, 0, alpha);
                        glVertex3f(positionX + h, positionY + h, positionZ + h);

                        //glColor4f(density011, 0, 0, alpha);
                        glVertex3f(positionX, positionY + h, positionZ + h);

                        //glColor4f(density001, 0, 0, alpha);
                        glVertex3f(positionX, positionY, positionZ + h);

                        //glColor4f(density101, 0, 0, alpha);
                        glVertex3f(positionX + h, positionY, positionZ + h);

                        //glColor4f(density110, 0, 0, alpha);
                        glVertex3f(positionX + h, positionY + h, positionZ);

                        //glColor4f(density111, 0, 0, alpha);
                        glVertex3f(positionX + h, positionY + h, positionZ + h);

                        //glColor4f(density101, 0, 0, alpha);
                        glVertex3f(positionX + h, positionY, positionZ + h);

                        //glColor4f(density100, 0, 0, alpha);
                        glVertex3f(positionX + h, positionY, positionZ);

                        //glColor4f(density010, 0, 0, alpha);
                        glVertex3f(positionX, positionY + h, positionZ);

                        //glColor4f(density110, 0, 0, alpha);
                        glVertex3f(positionX + h, positionY + h, positionZ);

                        //glColor4f(density100, 0, 0, alpha);
                        glVertex3f(positionX + h, positionY, positionZ);

                        //glColor4f(density000, 0, 0, alpha);
                        glVertex3f(positionX, positionY, positionZ);

                        //glColor4f(density011, 0, 0, alpha);
                        glVertex3f(positionX, positionY + h, positionZ + h);

                        //glColor4f(density010, 0, 0, alpha);
                        glVertex3f(positionX, positionY + h, positionZ);

                        //glColor4f(density000, 0, 0, alpha);
                        glVertex3f(positionX, positionY, positionZ);

                        //glColor4f(density001, 0, 0, alpha);
                        glVertex3f(positionX, positionY, positionZ + h);

                        //glColor4f(density100, 0, 0, alpha);
                        glVertex3f(positionX + h, positionY, positionZ);

                        //glColor4f(density000, 0, 0, alpha);
                        glVertex3f(positionX, positionY, positionZ);

                        //glColor4f(density001, 0, 0, alpha);
                        glVertex3f(positionX, positionY, positionZ + h);

                        //glColor4f(density101, 0, 0, alpha);
                        glVertex3f(positionX + h, positionY, positionZ + h);

                        //glColor4f(density110, 0, 0, alpha);
                        glVertex3f(positionX + h, positionY + h, positionZ);

                        //glColor4f(density010, 0, 0, alpha);
                        glVertex3f(positionX, positionY + h, positionZ);

                        //glColor4f(density011, 0, 0, alpha);
                        glVertex3f(positionX, positionY + h, positionZ + h);

                        //glColor4f(density111, 0, 0, alpha);
                        glVertex3f(positionX + h, positionY + h, positionZ + h);
                    }
                }

            }
        }
    }
    glEnd();
}

void Visualizer::wuDrawDensity()
{
    GLfloat positionX;
    GLfloat positionY;
    GLfloat positionZ;

    //GLfloat density000;
    //GLfloat density010;
    //GLfloat density100;
    //GLfloat density110;
    //GLfloat density001;
    //GLfloat density011;
    //GLfloat density101;
    //GLfloat density111;

    const GLfloat scale = 1.3f;
    const GLfloat h = scale / maxN;

    glBegin(GL_QUADS);

    for (int x = 1; x <= nx - 1; x++)
    {
        positionX = (x - 0.5f) * h - scale / 2.f;

        for (int y = 1; y <= ny - 1; y++)
        {
            positionY = (y - 0.5f) * h;

            for (int z = 1; z <= nz - 1; z++)
            {
                positionZ = (z - 0.5f) * h - scale / 2.f;
                
                //if (!solver->getObstacle(x,y,z))
                {
                    float dsB[2][2][2] = { 0.f, };
                    float dsR[2][2][2] = { 0.f, };

                    for (int ii = 0; ii < 2; ii++)
                        for (int jj = 0; jj < 2; jj++)
                            for (int kk = 0; kk < 2; kk++)
                            {
                                //float t = solver->getDensityGPU(x + ii, y + jj, z + kk);
                                float t = (solver->getDensityGPU(x + ii, y + jj, z + kk) - setting->eulerianSetting.minTem)
                                    / (setting->eulerianSetting.maxTem - setting->eulerianSetting.minTem);
                                //printf("t : %f \n", t);
                                dsR[ii][jj][jj] = t;
                                dsB[ii][jj][jj] = 1 - t;
                            }

                    glColor4f(dsR[1][1][1], 0, dsB[1][1][1], alpha);
                    glVertex3f(positionX + h, positionY + h, positionZ + h);

                    glColor4f(dsR[0][1][1], 0, dsB[0][1][1], alpha);
                    glVertex3f(positionX, positionY + h, positionZ + h);

                    glColor4f(dsR[0][0][1], 0, dsB[0][0][1], alpha);
                    glVertex3f(positionX, positionY, positionZ + h);

                    glColor4f(dsR[1][0][1], 0, dsB[1][0][1], alpha);
                    glVertex3f(positionX + h, positionY, positionZ + h);


                    glColor4f(dsR[1][1][0], 0, dsB[1][1][0], alpha);
                    glVertex3f(positionX + h, positionY + h, positionZ);

                    glColor4f(dsR[1][1][1], 0, dsB[1][1][1], alpha);
                    glVertex3f(positionX + h, positionY + h, positionZ + h);

                    glColor4f(dsR[1][0][1], 0, dsB[1][0][1], alpha);
                    glVertex3f(positionX + h, positionY, positionZ + h);

                    glColor4f(dsR[1][0][0], 0, dsB[1][0][0], alpha);
                    glVertex3f(positionX + h, positionY, positionZ);


                    glColor4f(dsR[0][1][0], 0, dsB[0][1][0], alpha);
                    glVertex3f(positionX, positionY + h, positionZ);

                    glColor4f(dsR[1][1][0], 0, dsB[1][1][0], alpha);
                    glVertex3f(positionX + h, positionY + h, positionZ);

                    glColor4f(dsR[1][0][0], 0, dsB[1][0][0], alpha);
                    glVertex3f(positionX + h, positionY, positionZ);

                    glColor4f(dsR[0][0][0], 0, dsB[0][0][0], alpha);
                    glVertex3f(positionX, positionY, positionZ);


                    glColor4f(dsR[0][1][1], 0, dsB[0][1][1], alpha);
                    glVertex3f(positionX, positionY + h, positionZ + h);

                    glColor4f(dsR[0][1][0], 0, dsB[0][1][0], alpha);
                    glVertex3f(positionX, positionY + h, positionZ);

                    glColor4f(dsR[0][0][0], 0, dsB[0][0][0], alpha);
                    glVertex3f(positionX, positionY, positionZ);

                    glColor4f(dsR[0][0][1], 0, dsB[0][0][1], alpha);
                    glVertex3f(positionX, positionY, positionZ + h);


                    glColor4f(dsR[1][0][0], 0, dsB[1][0][0], alpha);
                    glVertex3f(positionX + h, positionY, positionZ);

                    glColor4f(dsR[0][0][0], 0, dsB[0][0][0], alpha);
                    glVertex3f(positionX, positionY, positionZ);

                    glColor4f(dsR[0][0][1], 0, dsB[0][0][1], alpha);
                    glVertex3f(positionX, positionY, positionZ + h);

                    glColor4f(dsR[1][0][1], 0, dsB[1][0][1], alpha);
                    glVertex3f(positionX + h, positionY, positionZ + h);


                    glColor4f(dsR[1][1][0], 0, dsB[1][1][0], alpha);
                    glVertex3f(positionX + h, positionY + h, positionZ);

                    glColor4f(dsR[0][1][0], 0, dsB[0][1][0], alpha);
                    glVertex3f(positionX, positionY + h, positionZ);

                    glColor4f(dsR[0][1][1], 0, dsB[0][1][1], alpha);
                    glVertex3f(positionX, positionY + h, positionZ + h);

                    glColor4f(dsR[1][1][1], 0, dsB[1][1][1], alpha);
                    glVertex3f(positionX + h, positionY + h, positionZ + h);
                }
            }
        }
    }
    glEnd();
}

void Visualizer::wuTemptoGLColor(float temp, float opacity)
{
    opacity *= 0.7f;
    float t = (temp - minTem) / (maxTem - minTem);
    
    //if(false)
    {
        //if (opacity >= 1.5)
        //    glColor4f(1.f, 0.f, 0.f, 1.f);
        //else
        glColor4f(1.f, 1.f, 1.f, opacity);
        return;
    }

    if (t < 0) t = 0;
    if (t > 1) t = 1;

    if (t < 0.4f) {
        glColor4f(0, 0, 1.f, opacity);
    }
    else if (t < 0.6f) {
        glColor4f(0, 1 / 0.2f * (t - 0.4f), 1.f, opacity);
    }
    else if (t < 0.75f) {
        glColor4f(0, 1.f, 1 - 1 / 0.15f * (t - 0.6f), opacity);
    }
    else if (t < 0.8f) {
        glColor4f(1.f / 0.05f * (t - 0.75f), 1.f, 0, opacity);
    }
    else if (t < 0.85f) {
        glColor4f(1.f, 1 - 1 / 0.05f * (t - 0.8f), 0, opacity);
    }
    else {
        glColor4f(1.f, 0, 0, opacity);
    }
}

void Visualizer::wuDrawVelocity()
{
    GLfloat density000;
    GLfloat density010;
    GLfloat density100;
    GLfloat density110;
    GLfloat density001;
    GLfloat density011;
    GLfloat density101;
    GLfloat density111;

    GLfloat positionX;
    GLfloat positionY;
    GLfloat positionZ;
    GLfloat scale = 1.3f;
    GLfloat h = scale / maxN;

    glColor3f(1.f, 1.f, 1.f);



    for (int x = 1; x <= nx; x++)
    {
        positionX = (x - 0.5f) * h - scale / 2.f;
        for (int y = 1; y <= ny; y++)
        {
            positionY = (y - 0.5f) * h;
            for (int z = 1; z <= nz; z++)
            {
                positionZ = (z - 0.5f) * h - scale / 2.f;
                
                auto vu = solver->getVelocityGPUU(x, y, z) * drawScale;
                auto vv = solver->getVelocityGPUV(x, y, z) * drawScale;
                auto vw = solver->getVelocityGPUW(x, y, z) * drawScale;
                auto length = sqrt(vu * vu + vv * vv + vw * vw);


                if (length < 0.01f)
                    continue;

                if (vv < 0)
                    glColor3f(1.f, 1.f, 1.f);
                else
                    glColor3f(0.f, 1.0f, 0.f);

                glBegin(GL_LINES);
                glVertex3f(positionX - vu / 2, positionY - vv / 2, positionZ - vw / 2);
                glVertex3f(positionX + vu / 2, positionY + vv / 2, positionZ + vw / 2);
                glEnd();
            }
        }
    }

}

float Visualizer::wuGetDensity(float px, float py, float pz) {

    wu_no_stick(px, py, pz);

    int x, y, z;
    float xr, yr, zr;
    x = floor(px - 0.5f);
    y = floor(py - 0.5f);
    z = floor(pz - 0.5f);
    xr = px - 0.5f - x;
    yr = py - 0.5f - y;
    zr = pz - 0.5f - z;

    if (px < 0.5f) {
        x = 0;
        xr = 0;
    }
    if (py < 0.5f) {
        y = 0;
        yr = 0;
    }
    if (pz < 0.5f) {
        z = 0;
        zr = 0;
    }
    if (px >= nx + 1.5) {
        x = nx;
        xr = 1;
    }
    if (py >= ny + 1.5) {
        y = ny;
        yr = 1;
    }
    if (pz >= nz + 1.5) {
        z = nz;
        zr = 1;
    }

    float shape[8];
    for (int j = 0; j < 8; j++) {
        shape[j] = 1;
        shape[j] *= j % 2 ? xr : (1 - xr);
        shape[j] *= (j / 2 % 2) ? yr : (1 - yr);
        shape[j] *= (j / 4 % 2) ? zr : (1 - zr);
    }

    float d = 0;
    for (int j = 0; j < 8; j++) {
        d += shape[j] * solver->getDensityGPU(x + j % 2, y + (j / 2 % 2), z + (j / 4 % 2));
    }
    return d;
}

float Visualizer::wuVorticity(float px, float py, float pz) {

    int x, y, z;
    float xr, yr, zr;
    x = floor(px - 0.5f);
    y = floor(py - 0.5f);
    z = floor(pz - 0.5f);
    xr = px - 0.5f - x;
    yr = py - 0.5f - y;
    zr = pz - 0.5f - z;

    if (px < 0.5f) {
        x = 0;
        xr = 0;
    }
    if (py < 0.5f) {
        y = 0;
        yr = 0;
    }
    if (pz < 0.5f) {
        z = 0;
        zr = 0;
    }
    if (px >= nx - 1.5) {
        x = nx - 2;
        xr = 1;
    }
    if (py >= ny - 1.5) {
        y = ny - 2;
        yr = 1;
    }
    if (pz >= nz - 1.5) {
        z = nz - 2;
        zr = 1;
    }

    //vorticity
    
    float omega[3];
    {
        omega[0] = 0;

        omega[0] += xr * zr * (solver->getVelocityGPUW(x + 1, y + 1, z + 1)- solver->getVelocityGPUW(x + 1, y, z + 1));
        omega[0] += xr * (1 - zr) * (solver->getVelocityGPUW(x + 1, y + 1, z) - solver->getVelocityGPUW(x + 1, y, z));
        omega[0] += (1 - xr) * zr * (solver->getVelocityGPUW(x, y + 1, z + 1) - solver->getVelocityGPUW(x, y, z + 1));
        omega[0] += (1 - xr) * (1 - zr) * (solver->getVelocityGPUW(x, y + 1, z) - solver->getVelocityGPUW(x, y, z));

        omega[0] -= xr * yr * (solver->getVelocityGPUV(x + 1, y + 1, z + 1) - solver->getVelocityGPUV(x + 1, y + 1, z));
        omega[0] -= xr * (1 - yr) * (solver->getVelocityGPUV(x + 1, y, z + 1) - solver->getVelocityGPUV(x + 1, y, z));
        omega[0] -= (1 - xr) * yr * (solver->getVelocityGPUV(x, y + 1, z + 1) - solver->getVelocityGPUV(x, y + 1, z));
        omega[0] -= (1 - xr) * (1 - yr) * (solver->getVelocityGPUV(x, y, z + 1) - solver->getVelocityGPUV(x, y, z));
    }
    {
        omega[1] = 0;

        omega[1] += xr * yr * (solver->getVelocityGPUU(x + 1, y + 1, z + 1) - solver->getVelocityGPUU(x + 1, y + 1, z));
        omega[1] += xr * (1 - yr) * (solver->getVelocityGPUU(x + 1, y, z + 1) - solver->getVelocityGPUU(x + 1, y, z));
        omega[1] += (1 - xr) * yr * (solver->getVelocityGPUU(x, y + 1, z + 1) - solver->getVelocityGPUU(x, y + 1, z));
        omega[1] += (1 - xr) * (1 - yr) * (solver->getVelocityGPUU(x, y, z + 1) - solver->getVelocityGPUU(x, y, z));

        omega[1] -= yr * zr * (solver->getVelocityGPUW(x + 1, y + 1, z + 1) - solver->getVelocityGPUW(x, y + 1, z + 1));
        omega[1] -= yr * (1 - zr) * (solver->getVelocityGPUW(x + 1, y + 1, z) - solver->getVelocityGPUW(x, y + 1, z));
        omega[1] -= (1 - yr) * zr * (solver->getVelocityGPUW(x + 1, y, z + 1) - solver->getVelocityGPUW(x, y, z + 1));
        omega[1] -= (1 - yr) * (1 - zr) * (solver->getVelocityGPUW(x + 1, y, z) - solver->getVelocityGPUW(x, y, z));
    }
    {
        omega[2] = 0;

        omega[2] -= yr * zr * (solver->getVelocityGPUV(x + 1, y + 1, z + 1) - solver->getVelocityGPUV(x, y + 1, z + 1));
        omega[2] -= yr * (1 - zr) * (solver->getVelocityGPUV(x + 1, y + 1, z) - solver->getVelocityGPUV(x, y + 1, z));
        omega[2] -= (1 - yr) * zr * (solver->getVelocityGPUV(x + 1, y, z + 1) - solver->getVelocityGPUV(x, y, z + 1));
        omega[2] -= (1 - yr) * (1 - zr) * (solver->getVelocityGPUV(x + 1, y, z) - solver->getVelocityGPUV(x, y, z));

        omega[2] += xr * zr * (solver->getVelocityGPUU(x + 1, y + 1, z + 1) - solver->getVelocityGPUU(x + 1, y, z + 1));
        omega[2] += xr * (1 - zr) * (solver->getVelocityGPUU(x + 1, y + 1, z) - solver->getVelocityGPUU(x + 1, y, z));
        omega[2] += (1 - xr) * zr * (solver->getVelocityGPUU(x, y + 1, z + 1) - solver->getVelocityGPUU(x, y, z + 1));
        omega[2] += (1 - xr) * (1 - zr) * (solver->getVelocityGPUU(x, y + 1, z) - solver->getVelocityGPUU(x, y, z));
    }
    //vorticity

    return sqrt(omega[0] * omega[0] + omega[1] * omega[1] + omega[2] * omega[2]);
}

bool Visualizer::wuInVortex(float px, float py, float pz) {

    int x, y, z;
    float xr, yr, zr;
    x = floor(px - 0.5f);
    y = floor(py - 0.5f);
    z = floor(pz - 0.5f);
    xr = px - 0.5f - x;
    yr = py - 0.5f - y;
    zr = pz - 0.5f - z;

    if (px < 0.5f) {
        x = 0;
        xr = 0;
    }
    if (py < 0.5f) {
        y = 0;
        yr = 0;
    }
    if (pz < 0.5f) {
        z = 0;
        zr = 0;
    }
    if (px >= nx - 1.5) {
        x = nx - 2;
        xr = 1;
    }
    if (py >= ny - 1.5) {
        y = ny - 2;
        yr = 1;
    }
    if (pz >= nz - 1.5) {
        z = nz - 2;
        zr = 1;
    }

    //vorticity

    float J[3][3];
    {
        J[0][0] = 0;

        J[0][0] += yr * zr * (solver->getVelocityGPUU(x + 1, y + 1, z + 1) - solver->getVelocityGPUU(x, y + 1, z + 1));
        J[0][0] += yr * (1 - zr) * (solver->getVelocityGPUU(x + 1, y + 1, z) - solver->getVelocityGPUU(x, y + 1, z));
        J[0][0] += (1 - yr) * zr * (solver->getVelocityGPUU(x + 1, y, z + 1) - solver->getVelocityGPUU(x, y, z + 1));
        J[0][0] += (1 - yr) * (1 - zr) * (solver->getVelocityGPUU(x + 1, y, z) - solver->getVelocityGPUU(x, y, z));

        J[1][0] = 0;

        J[1][0] += yr * zr * (solver->getVelocityGPUV(x + 1, y + 1, z + 1) - solver->getVelocityGPUV(x, y + 1, z + 1));
        J[1][0] += yr * (1 - zr) * (solver->getVelocityGPUV(x + 1, y + 1, z) - solver->getVelocityGPUV(x, y + 1, z));
        J[1][0] += (1 - yr) * zr * (solver->getVelocityGPUV(x + 1, y, z + 1) - solver->getVelocityGPUV(x, y, z + 1));
        J[1][0] += (1 - yr) * (1 - zr) * (solver->getVelocityGPUV(x + 1, y, z) - solver->getVelocityGPUV(x, y, z));

        J[2][0] = 0;

        J[2][0] += yr * zr * (solver->getVelocityGPUW(x + 1, y + 1, z + 1) - solver->getVelocityGPUW(x, y + 1, z + 1));
        J[2][0] += yr * (1 - zr) * (solver->getVelocityGPUW(x + 1, y + 1, z) - solver->getVelocityGPUW(x, y + 1, z));
        J[2][0] += (1 - yr) * zr * (solver->getVelocityGPUW(x + 1, y, z + 1) - solver->getVelocityGPUW(x, y, z + 1));
        J[2][0] += (1 - yr) * (1 - zr) * (solver->getVelocityGPUW(x + 1, y, z) - solver->getVelocityGPUW(x, y, z));
    }
    {
        J[0][0] = 0;

        J[0][0] += yr * zr * (solver->getVelocityGPUU(x + 1, y + 1, z + 1) - solver->getVelocityGPUU(x, y + 1, z + 1));
        J[0][0] += yr * (1 - zr) * (solver->getVelocityGPUU(x + 1, y + 1, z) - solver->getVelocityGPUU(x, y + 1, z));
        J[0][0] += (1 - yr) * zr * (solver->getVelocityGPUU(x + 1, y, z + 1) - solver->getVelocityGPUU(x, y, z + 1));
        J[0][0] += (1 - yr) * (1 - zr) * (solver->getVelocityGPUU(x + 1, y, z) - solver->getVelocityGPUU(x, y, z));

        J[1][0] = 0;

        J[1][0] += yr * zr * (solver->getVelocityGPUV(x + 1, y + 1, z + 1) - solver->getVelocityGPUV(x, y + 1, z + 1));
        J[1][0] += yr * (1 - zr) * (solver->getVelocityGPUV(x + 1, y + 1, z) - solver->getVelocityGPUV(x, y + 1, z));
        J[1][0] += (1 - yr) * zr * (solver->getVelocityGPUV(x + 1, y, z + 1) - solver->getVelocityGPUV(x, y, z + 1));
        J[1][0] += (1 - yr) * (1 - zr) * (solver->getVelocityGPUV(x + 1, y, z) - solver->getVelocityGPUV(x, y, z));

        J[2][0] = 0;

        J[2][0] += yr * zr * (solver->getVelocityGPUW(x + 1, y + 1, z + 1) - solver->getVelocityGPUW(x, y + 1, z + 1));
        J[2][0] += yr * (1 - zr) * (solver->getVelocityGPUW(x + 1, y + 1, z) - solver->getVelocityGPUW(x, y + 1, z));
        J[2][0] += (1 - yr) * zr * (solver->getVelocityGPUW(x + 1, y, z + 1) - solver->getVelocityGPUW(x, y, z + 1));
        J[2][0] += (1 - yr) * (1 - zr) * (solver->getVelocityGPUW(x + 1, y, z) - solver->getVelocityGPUW(x, y, z));
    }
    {
        J[0][1] = 0;

        J[0][1] += xr * zr * (solver->getVelocityGPUU(x + 1, y + 1, z + 1) - solver->getVelocityGPUU(x + 1, y, z + 1));
        J[0][1] += xr * (1 - zr) * (solver->getVelocityGPUU(x + 1, y + 1, z) - solver->getVelocityGPUU(x + 1, y, z));
        J[0][1] += (1 - xr) * zr * (solver->getVelocityGPUU(x, y + 1, z + 1) - solver->getVelocityGPUU(x, y, z + 1));
        J[0][1] += (1 - xr) * (1 - zr) * (solver->getVelocityGPUU(x, y + 1, z) - solver->getVelocityGPUU(x, y, z));

        J[1][1] = 0;

        J[1][1] += xr * zr * (solver->getVelocityGPUV(x + 1, y + 1, z + 1) - solver->getVelocityGPUV(x + 1, y, z + 1));
        J[1][1] += xr * (1 - zr) * (solver->getVelocityGPUV(x + 1, y + 1, z) - solver->getVelocityGPUV(x + 1, y, z));
        J[1][1] += (1 - xr) * zr * (solver->getVelocityGPUV(x, y + 1, z + 1) - solver->getVelocityGPUV(x, y, z + 1));
        J[1][1] += (1 - xr) * (1 - zr) * (solver->getVelocityGPUV(x, y + 1, z) - solver->getVelocityGPUV(x, y, z));

        J[2][1] = 0;

        J[2][1] += xr * zr * (solver->getVelocityGPUW(x + 1, y + 1, z + 1) - solver->getVelocityGPUW(x + 1, y, z + 1));
        J[2][1] += xr * (1 - zr) * (solver->getVelocityGPUW(x + 1, y + 1, z) - solver->getVelocityGPUW(x + 1, y, z));
        J[2][1] += (1 - xr) * zr * (solver->getVelocityGPUW(x, y + 1, z + 1) - solver->getVelocityGPUW(x, y, z + 1));
        J[2][1] += (1 - xr) * (1 - zr) * (solver->getVelocityGPUW(x, y + 1, z) - solver->getVelocityGPUW(x, y, z));
    }

    {
        J[0][2] = 0;

        J[0][2] += xr * yr * (solver->getVelocityGPUU(x + 1, y + 1, z + 1) - solver->getVelocityGPUU(x + 1, y + 1, z));
        J[0][2] += xr * (1 - yr) * (solver->getVelocityGPUU(x + 1, y, z + 1) - solver->getVelocityGPUU(x + 1, y, z));
        J[0][2] += (1 - xr) * yr * (solver->getVelocityGPUU(x, y + 1, z + 1) - solver->getVelocityGPUU(x, y + 1, z));
        J[0][2] += (1 - xr) * (1 - yr) * (solver->getVelocityGPUU(x, y, z + 1) - solver->getVelocityGPUU(x, y, z));

        J[1][2] = 0;

        J[1][2] += xr * yr * (solver->getVelocityGPUV(x + 1, y + 1, z + 1) - solver->getVelocityGPUV(x + 1, y + 1, z));
        J[1][2] += xr * (1 - yr) * (solver->getVelocityGPUV(x + 1, y, z + 1) - solver->getVelocityGPUV(x + 1, y, z));
        J[1][2] += (1 - xr) * yr * (solver->getVelocityGPUV(x, y + 1, z + 1) - solver->getVelocityGPUV(x, y + 1, z));
        J[1][2] += (1 - xr) * (1 - yr) * (solver->getVelocityGPUV(x, y, z + 1) - solver->getVelocityGPUV(x, y, z));

        J[2][2] = 0;

        J[2][2] += xr * yr * (solver->getVelocityGPUW(x + 1, y + 1, z + 1) - solver->getVelocityGPUW(x + 1, y + 1, z));
        J[2][2] += xr * (1 - yr) * (solver->getVelocityGPUW(x + 1, y, z + 1) - solver->getVelocityGPUW(x + 1, y, z));
        J[2][2] += (1 - xr) * yr * (solver->getVelocityGPUW(x, y + 1, z + 1) - solver->getVelocityGPUW(x, y + 1, z));
        J[2][2] += (1 - xr) * (1 - yr) * (solver->getVelocityGPUW(x, y, z + 1) - solver->getVelocityGPUW(x, y, z));
    }

    float S[3][3];
    float O[3][3];

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            S[i][j] = (J[i][j] + J[j][i]) / 2;
            O[i][j] = (J[i][j] - J[j][i]) / 2;
        }

    float Arr[3][3];

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            Arr[i][j] = S[i][0] * S[0][j] + S[i][1] * S[1][j] + S[i][2] * S[2][j]\
                + O[i][0] * O[0][j] + O[i][1] * O[1][j] + O[i][2] * O[2][j];
        }

    //cubic solver
    float b, c, d;
    b = -(Arr[0][0] + Arr[1][1] + Arr[2][2]);
    c = (Arr[0][0] * (Arr[1][1] + Arr[2][2]) + Arr[1][1] * Arr[2][2]) - (Arr[0][1] * Arr[0][1] + Arr[0][2] * Arr[0][2] + Arr[1][2] * Arr[1][2]);
    d = -(Arr[0][0] * Arr[1][1] * Arr[2][2] + 2 * Arr[0][1] * Arr[0][2] * Arr[1][2]) + (Arr[0][0] * Arr[1][2] * Arr[1][2] + Arr[1][1] * Arr[0][2] * Arr[0][2] + Arr[2][2] * Arr[0][1] * Arr[0][1]);

    float q, disc, r, s, t, term1;
    term1 = b / 3;
    q = (3 * c - b * b) / 9;
    r = (-27 * d + b * (9 * c - 2 * b * b)) / 54;
    disc = q * q * q + r * r;
    if (disc > 0) {
        printf("Complex Eigenvalues\n");
        return false;
    }
    //s = r + sqrt(disc);
    //t = r - sqrt(disc);
    //s = (s < 0) ? (-pow(-s, 1.0 / 3)) : (pow(s, 1.0 / 3));
    //t = (s < 0) ? (-pow(-s, 1.0 / 3)) : (pow(s, 1.0 / 3));

    float dum1, r13;
    q = -q;
    dum1 = q * q * q;
    dum1 = acos(r / sqrt(dum1));
    r13 = 2 * sqrt(q);


    float l1, l2, l3;
    l1 = -term1 + r13 * cos(dum1 / 3);
    l2 = -term1 + r13 * cos((dum1 + 2 * 3.14159265358979) / 3);
    l3 = -term1 + r13 * cos((dum1 + 4 * 3.14159265358979) / 3);

    int abovz = 0;
    if (l1 > 0) abovz++;
    if (l2 > 0) abovz++;
    if (l3 > 0) abovz++;

    if (abovz < 2) return true;

    return false;
}


void Visualizer::wuGetVelocity(float px, float py, float pz, float& vx, float& vy, float& vz) {

    int x, y, z;
    float xr, yr, zr;
    x = floor(px - 0.5f);
    y = floor(py - 0.5f);
    z = floor(pz - 0.5f);
    xr = px - 0.5f - x;
    yr = py - 0.5f - y;
    zr = pz - 0.5f - z;

    if (px < 0.5f) {
        x = 0;
        xr = 0;
    }
    if (py < 0.5f) {
        y = 0;
        yr = 0;
    }
    if (pz < 0.5f) {
        z = 0;
        zr = 0;
    }
    if (px >= nx + 1.5) {
        x = nx;
        xr = 1;
    }
    if (py >= ny + 1.5) {
        y = ny;
        yr = 1;
    }
    if (pz >= nz + 1.5) {
        z = nz;
        zr = 1;
    }
    
    float shape[8];
    for (int j = 0; j < 8; j++) {
        shape[j] = 1;
        shape[j] *= j % 2 ? xr : (1 - xr);
        shape[j] *= (j / 2 % 2) ? yr : (1 - yr);
        shape[j] *= (j / 4 % 2) ? zr : (1 - zr);
    }

    vx = vy = vz = 0;
    for (int j = 0; j < 8; j++) {
        vx += shape[j] * solver->getVelocityGPUU(x + j % 2, y + (j / 2 % 2), z + (j / 4 % 2));
        vy += shape[j] * solver->getVelocityGPUV(x + j % 2, y + (j / 2 % 2), z + (j / 4 % 2));
        vz += shape[j] * solver->getVelocityGPUW(x + j % 2, y + (j / 2 % 2), z + (j / 4 % 2));
    }
}

void Visualizer::wuGetVelocityRK4(float px, float py, float pz, float& vx, float& vy, float& vz)
{
    float dp[4][3];
    /*
    float x, y, z;
    x = px; y = py; z = pz;
    wuGetVelocity(x, y, z, dp[0][0], dp[0][1], dp[0][2]);

    x = px + particleSpeed * dp[0][0] * 0.5f;
    y = py + particleSpeed * dp[0][1] * 0.5f;
    z = pz + particleSpeed * dp[0][2] * 0.5f;
    wu_no_stick(x, y, z);
    wuGetVelocity(x, y, z, dp[1][0], dp[1][1], dp[1][2]);

    x = px + particleSpeed * dp[1][0] * 0.5f;
    y = py + particleSpeed * dp[1][1] * 0.5f;
    z = pz + particleSpeed * dp[1][2] * 0.5f;
    wu_no_stick(x, y, z);
    wuGetVelocity(x, y, z, dp[2][0], dp[2][1], dp[2][2]);

    x = px + particleSpeed * dp[2][0];
    y = py + particleSpeed * dp[2][1];
    z = pz + particleSpeed * dp[2][2];
    wu_no_stick(x, y, z);
    wuGetVelocity(x, y, z, dp[3][0], dp[3][1], dp[3][2]);
    */

    wuGetVelocity(px, py, pz, dp[0][0], dp[0][1], dp[0][2]);
    wuGetVelocity(px + particleSpeed * dp[0][0] * 0.5f, py + particleSpeed * dp[0][1] * 0.5f, pz + particleSpeed * dp[0][2] * 0.5f, dp[1][0], dp[1][1], dp[1][2]);
    wuGetVelocity(px + particleSpeed * dp[1][0] * 0.5f, py + particleSpeed * dp[1][1] * 0.5f, pz + particleSpeed * dp[1][2] * 0.5f, dp[2][0], dp[2][1], dp[2][2]);
    wuGetVelocity(px + particleSpeed * dp[2][0], py + particleSpeed * dp[2][1], pz + particleSpeed * dp[2][2], dp[3][0], dp[3][1], dp[3][2]);



    vx = dp[0][0] / 6.f + dp[1][0] / 3.f + dp[2][0] / 3.f + dp[3][0] / 6.f;
    vy = dp[0][1] / 6.f + dp[1][1] / 3.f + dp[2][1] / 3.f + dp[3][1] / 6.f;
    vz = dp[0][2] / 6.f + dp[1][2] / 3.f + dp[2][2] / 3.f + dp[3][2] / 6.f;

}

void Visualizer::wuGetNextStreamPosition(float px, float py, float pz, float& vx, float& vy, float& vz)
{
    int ix = floor(px);
    int iy = floor(py);
    int iz = floor(pz);
    float rx = px - ix;
    float ry = py - iy;
    float rz = pz - iz;

    float W0[3], W1[3], W2[3], W3[3];
    float V[3][3];

    W0[0] = solver->getVelocityGPUU(ix, iy, iz);
    W0[1] = solver->getVelocityGPUV(ix, iy, iz);
    W0[2] = solver->getVelocityGPUW(ix, iy, iz);
    W3[0] = solver->getVelocityGPUU(ix + 1, iy + 1, iz + 1);
    W3[1] = solver->getVelocityGPUV(ix + 1, iy + 1, iz + 1);
    W3[2] = solver->getVelocityGPUW(ix + 1, iy + 1, iz + 1);

    if (rx >= ry && ry >= rz) {
        W1[0] = solver->getVelocityGPUU(ix + 1, iy, iz);
        W1[1] = solver->getVelocityGPUV(ix + 1, iy, iz);
        W1[2] = solver->getVelocityGPUW(ix + 1, iy, iz);
        W2[0] = solver->getVelocityGPUU(ix + 1, iy + 1, iz);
        W2[1] = solver->getVelocityGPUV(ix + 1, iy + 1, iz);
        W2[2] = solver->getVelocityGPUW(ix + 1, iy + 1, iz);
        for (int i = 0; i < 3; i++) {
            V[i][0] = W1[i] - W0[i];
            V[i][1] = W2[i] - W1[i];
            V[i][2] = W3[i] - W2[i];
        }
    }
    else if (rx > rz && rz >= ry) {
        W1[0] = solver->getVelocityGPUU(ix + 1, iy, iz);
        W1[1] = solver->getVelocityGPUV(ix + 1, iy, iz);
        W1[2] = solver->getVelocityGPUW(ix + 1, iy, iz);
        W2[0] = solver->getVelocityGPUU(ix + 1, iy, iz + 1);
        W2[1] = solver->getVelocityGPUV(ix + 1, iy, iz + 1);
        W2[2] = solver->getVelocityGPUW(ix + 1, iy, iz + 1);
        for (int i = 0; i < 3; i++) {
            V[i][0] = W1[i] - W0[i];
            V[i][2] = W2[i] - W1[i];
            V[i][1] = W3[i] - W2[i];
        }
    }
    else if (ry >= rx && rx > rz) {
        W1[0] = solver->getVelocityGPUU(ix, iy + 1, iz);
        W1[1] = solver->getVelocityGPUV(ix, iy + 1, iz);
        W1[2] = solver->getVelocityGPUW(ix, iy + 1, iz);
        W2[0] = solver->getVelocityGPUU(ix + 1, iy + 1, iz);
        W2[1] = solver->getVelocityGPUV(ix + 1, iy + 1, iz);
        W2[2] = solver->getVelocityGPUW(ix + 1, iy + 1, iz);
        for (int i = 0; i < 3; i++) {
            V[i][1] = W1[i] - W0[i];
            V[i][0] = W2[i] - W1[i];
            V[i][2] = W3[i] - W2[i];
        }
    }
    else if (ry > rz && rz >= rx) {
        W1[0] = solver->getVelocityGPUU(ix, iy + 1, iz);
        W1[1] = solver->getVelocityGPUV(ix, iy + 1, iz);
        W1[2] = solver->getVelocityGPUW(ix, iy + 1, iz);
        W2[0] = solver->getVelocityGPUU(ix, iy + 1, iz + 1);
        W2[1] = solver->getVelocityGPUV(ix, iy + 1, iz + 1);
        W2[2] = solver->getVelocityGPUW(ix, iy + 1, iz + 1);
        for (int i = 0; i < 3; i++) {
            V[i][1] = W1[i] - W0[i];
            V[i][2] = W2[i] - W1[i];
            V[i][0] = W3[i] - W2[i];
        }
    }
    else if (rz >= rx && rx > ry) {
        W1[0] = solver->getVelocityGPUU(ix, iy, iz + 1);
        W1[1] = solver->getVelocityGPUV(ix, iy, iz + 1);
        W1[2] = solver->getVelocityGPUW(ix, iy, iz + 1);
        W2[0] = solver->getVelocityGPUU(ix + 1, iy, iz + 1);
        W2[1] = solver->getVelocityGPUV(ix + 1, iy, iz + 1);
        W2[2] = solver->getVelocityGPUW(ix + 1, iy, iz + 1);
        for (int i = 0; i < 3; i++) {
            V[i][2] = W1[i] - W0[i];
            V[i][0] = W2[i] - W1[i];
            V[i][1] = W3[i] - W2[i];
        }
    }
    else if (rz > ry && ry > rx) {
        W1[0] = solver->getVelocityGPUU(ix, iy, iz + 1);
        W1[1] = solver->getVelocityGPUV(ix, iy, iz + 1);
        W1[2] = solver->getVelocityGPUW(ix, iy, iz + 1);
        W2[0] = solver->getVelocityGPUU(ix, iy + 1, iz + 1);
        W2[1] = solver->getVelocityGPUV(ix, iy + 1, iz + 1);
        W2[2] = solver->getVelocityGPUW(ix, iy + 1, iz + 1);
        for (int i = 0; i < 3; i++) {
            V[i][2] = W1[i] - W0[i];
            V[i][1] = W2[i] - W1[i];
            V[i][0] = W3[i] - W2[i];
        }
    }

    
    float l0, l1, l2;
    float ev0[3], ev1[3], ev2[3];

    



}

bool Visualizer::wu_no_stick(float& x, float& y, float& z) {
    const float no_stick = 0.001;
    bool ret = false;

    if (x < 0) {
        x = 0 + no_stick;
        ret = true;
    }
    if (y < 0) {
        y = 0 + no_stick;
        ret = true;
    }
    if (z < 0) {
        z = 0 + no_stick;
        ret = true;
    }

    if (x > nx) {
        x = nx - no_stick;
        ret = true;
    }
    if (y > ny) {
        y = ny - no_stick;
        ret = true;
    }
    if (z > nz) {
        z = nz - no_stick;
        ret = true;
    }

    return ret;
}

void Visualizer::advectParticles()
{
    for (int i = 0; i < numParticles; i++) {
        int length = particles[i].traceLength < lineLength ? particles[i].traceLength : lineLength;
        if (length > 0)
            particles[i].age++;

        for (int j = 0; j < length; j++) {
            float dx, dy, dz;
            wuGetVelocity(particles[i].particlePos[j][0], particles[i].particlePos[j][1], particles[i].particlePos[j][2], dx, dy, dz);

            particles[i].particlePos[j][0] += dx * particleSpeed;
            particles[i].particlePos[j][1] += dy * particleSpeed;
            particles[i].particlePos[j][2] += dz * particleSpeed;

            wu_no_stick(particles[i].particlePos[j][0], particles[i].particlePos[j][1], particles[i].particlePos[j][2]);
        }
        
    }
}


void Visualizer::advectParticlesRK4()
{
    for (int i = 0; i < numParticles; i++) {
        int length = particles[i].traceLength < lineLength ? particles[i].traceLength : lineLength;
        if (length > 0)
            particles[i].age++;

        for (int j = 0; j < length; j++) {
            float dx, dy, dz;
            wuGetVelocityRK4(particles[i].particlePos[j][0], particles[i].particlePos[j][1], particles[i].particlePos[j][2], dx, dy, dz);

            particles[i].particlePos[j][0] += dx * particleSpeed;
            particles[i].particlePos[j][1] += dy * particleSpeed;
            particles[i].particlePos[j][2] += dz * particleSpeed;

            wu_no_stick(particles[i].particlePos[j][0], particles[i].particlePos[j][1], particles[i].particlePos[j][2]);
        }

    }
}

int HashBig(int x, int n)
{
    long long A = x * 73856093;
    long long B = x * 19349663;
    long long C = x * 83492791;
    A = A ^ B ^ C;
    int r = (int)(A % n);
    if (r < 0) r += n;
    return r;
}


#define stLen 230
void Visualizer::wuDrawStreamline(bool animated)
{
    GLfloat positionX;
    GLfloat positionY;
    GLfloat positionZ;
    GLfloat prevPositionX;
    GLfloat prevPositionY;
    GLfloat prevPositionZ;

    GLfloat scale = 1.3f;
    GLfloat h = scale / maxN;

    AirConditioner& ac = solver->airconditioners[0];
    int x, y, z;
    x = ac.pos[0];
    y = ac.pos[1];
    z = ac.pos[2];
    
    float seeds[140][3];
    int numSeed = 0;
    
    if(ac.acType == ACTYPE::Ceiling)
    {
        //17*4=68 seeds
        for (float i = -2; i <= 2; i += 0.25f) {
            seeds[numSeed][0] = x - 4.2f + 0.5f;
            seeds[numSeed][1] = y - 1 + 0.5f;
            seeds[numSeed][2] = z + i + 0.5f;
            numSeed++;

            seeds[numSeed][0] = x + 4.2f + 0.5f;
            seeds[numSeed][1] = y - 1 + 0.5f;
            seeds[numSeed][2] = z + i + 0.5f;
            numSeed++;

            seeds[numSeed][0] = x + i + 0.5f;
            seeds[numSeed][1] = y - 1 + 0.5f;
            seeds[numSeed][2] = z - 4.2f + 0.5f;
            numSeed++;

            seeds[numSeed][0] = x + i + 0.5f;
            seeds[numSeed][1] = y - 1 + 0.5f;
            seeds[numSeed][2] = z + 4.2f + 0.5f;
            numSeed++;

        }
    }
    else if (ac.acType == ACTYPE::Tower) {
        //printf("Tower");
        if (ac.direction % 2 == 0) {
            //17*4+16=84 seeds

            int horX = 0, horZ = 0; //right-end direction
            float dx = 0, dz = 0;
            if (ac.direction == 0) {
                horX = 0;
                horZ = 1;
                dx = 0.2f;
                dz = 0;
            }
            else if (ac.direction == 2) {
                horX = 1;
                horZ = 0;
                dx = 0;
                dz = -0.2f;
            }
            else if (ac.direction == 4) {
                horX = 0;
                horZ = -1;
                dx = -0.2f;
                dz = 0;
            }
            else if (ac.direction == 6) {
                horX = -1;
                horZ = 0;
                dx = 0;
                dz = 0.2f;
            }

            for (float i = 0; i <= 4; i += 0.25f) {
                
                seeds[numSeed][0] = x - 2 * horX;
                seeds[numSeed][1] = y + 5 + i;
                seeds[numSeed][2] = z - 2 * horZ;
                numSeed++;

                seeds[numSeed][0] = x - 2 * horX;
                seeds[numSeed][1] = y + 10 + i;
                seeds[numSeed][2] = z - 2 * horZ;
                numSeed++;

                seeds[numSeed][0] = x + 1 * horX;
                seeds[numSeed][1] = y + 5 + i;
                seeds[numSeed][2] = z + 1 * horZ;
                numSeed++;

                seeds[numSeed][0] = x + 1 * horX;
                seeds[numSeed][1] = y + 10 + i;
                seeds[numSeed][2] = z + 1 * horZ;
                numSeed++;
            }

            for (int i = 0; i < numSeed; i++){
                seeds[i][0] += 0.5f + 0.2f * horZ;
                seeds[i][1] += 0.5f;
                seeds[i][2] += 0.5f - 0.2f * horX;
            }

            {
                float center[3];
                center[0] = x - 0.5f * horX + 0.5f;
                center[1] = y + 17.f;
                center[2] = z - 0.5f * horZ + 0.5f;
                float hor[3];
                hor[0] = -1.5f * horX;
                hor[1] = 0;
                hor[2] = -1.5f * horZ;

                int numRing = 16;
                for (int i = 0; i < numRing; i++) {
                    const float theta = i * 2 * 3.141592f / numRing;

                    seeds[numSeed][0] = center[0] + 1.0f * (hor[0] * cos(theta));
                    seeds[numSeed][1] = center[1] + 1.0f * 1.5f * sin(theta);
                    seeds[numSeed][2] = center[2] + 1.0f * (hor[2] * cos(theta));
                    numSeed++;
                }
            }
        }
        else {
            //17*4+16=84 seeds

            int horX = 0, horZ = 0; //right-end direction
            float dx = 0, dz = 0;
            if (ac.direction == 1) {
                horX = 1;
                horZ = 1;
                dx = 0.1f;
                dz = -0.1f;
            }
            else if (ac.direction == 3) {
                horX = 1;
                horZ = -1;
                dx = -0.1f;
                dz = -0.1f;
            }
            else if (ac.direction == 5) {
                horX = -1;
                horZ = -1;
                dx = -0.1f;
                dz = 0.1f;
            }
            else if (ac.direction == 7) {
                horX = -1;
                horZ = 1;
                dx = 0.1f;
                dz = 0.1f;
            }

            for (float i = 0; i <= 4; i += 0.25f) {

                seeds[numSeed][0] = x - horX;
                seeds[numSeed][1] = y + 5 + i;
                seeds[numSeed][2] = z - horZ;
                numSeed++;

                seeds[numSeed][0] = x - horX;
                seeds[numSeed][1] = y + 10 + i;
                seeds[numSeed][2] = z - horZ;
                numSeed++;

                seeds[numSeed][0] = x + horX;
                seeds[numSeed][1] = y + 5 + i;
                seeds[numSeed][2] = z + horZ;
                numSeed++;

                seeds[numSeed][0] = x + horX;
                seeds[numSeed][1] = y + 10 + i;
                seeds[numSeed][2] = z + horZ;
                numSeed++;
            }

            for (int i = 0; i < numSeed; i++) {
                seeds[i][0] += 0.5f + 0.1f * horZ;
                seeds[i][1] += 0.5f;
                seeds[i][2] += 0.5f - 0.1f * horX;
            }

            {
                float center[3];
                center[0] = x + 0.5f;
                center[1] = y + 17.f;
                center[2] = z + 0.5f;
                float hor[3];
                hor[0] = -1.5f / sqrt(2) * horX;
                hor[1] = 0;
                hor[2] = -1.5f / sqrt(2) * horZ;

                int numRing = 16;
                for (int i = 0; i < numRing; i++) {
                    const float theta = i * 2 * 3.141592f / numRing;

                    seeds[numSeed][0] = center[0] + 1.0f * (hor[0] * cos(theta));
                    seeds[numSeed][1] = center[1] + 1.0f * 1.5f * sin(theta);
                    seeds[numSeed][2] = center[2] + 1.0f * (hor[2] * cos(theta));
                    numSeed++;
                }
            }
        }
    }
    //printf("%d\n", numSeed);
    //for (int s = 0; s < numSeed; s++) {
    //    printf("(%f,%f,%f)\n", seeds[s][0], seeds[s][1], seeds[s][2]);
    //}

    bool isLineMode = true;

    glColor3f(1.f, 1.f, 1.f);
    if (isLineMode) {
        glBegin(GL_LINES);
    }
    else {
        glBegin(GL_POINTS);
    }
    for (int i = 0; i < numSeed; i++) {

        prevPositionX = seeds[i][0];
        prevPositionY = seeds[i][1];
        prevPositionZ = seeds[i][2];

        for (int j = 0; j < 180; j++) {

            float vx, vy, vz;
            wuGetVelocityRK4(prevPositionX, prevPositionY, prevPositionZ, vx, vy, vz);
            positionX = prevPositionX + particleSpeed * vx;
            positionY = prevPositionY + particleSpeed * vy;
            positionZ = prevPositionZ + particleSpeed * vz;
            wu_no_stick(positionX, positionY, positionZ);
            
            static const float lineSpeed = 50;
            static const int lineLength = 5;

            int turn = 5 * HashBig(i, 70) + j - (int)(solver->getFramestep() * 0.005f * lineSpeed);
            turn = turn % 70;
            if (turn < 0) turn += 70;

            if (turn < lineLength || !animated) {
                if (isLineMode) {
                    glVertex3f((positionX - 1) * h - scale / 2, (positionY - 1) * h, (positionZ - 1) * h - scale / 2);
                    glVertex3f((prevPositionX - 1) * h - scale / 2, (prevPositionY - 1) * h, (prevPositionZ - 1) * h - scale / 2);
                }
                else {
                    glVertex3f((positionX - 1)* h - scale / 2, (positionY - 1)* h, (positionZ - 1)* h - scale / 2);
                }

            }
            
            prevPositionX = positionX;
            prevPositionY = positionY;
            prevPositionZ = positionZ;

        }

    }
    glEnd();
}

float Dist(float* v0, float* v1) {
    return sqrt((v0[0] - v1[0]) * (v0[0] - v1[0]) + (v0[1] - v1[1]) * (v0[1] - v1[1]) + (v0[2] - v1[2]) * (v0[2] - v1[2]));
}
float Curv(float* u, float* v, float* w) {
    float d1 = Dist(u, v);
    float d2 = Dist(w, v);
    float d3 = 0;
    d3 += (u[0] - v[0]) * (w[0] - v[0]);
    d3 += (u[1] - v[1]) * (w[1] - v[1]);
    d3 += (u[2] - v[2]) * (w[2] - v[2]);
    return ((d3 / (d1 * d2)) + 1) / 2.0f;
}


float Visualizer::areaQuad(float* v0, float* v1, float* v2, float* v3) {

    float e1[3];
    float e2[3];
    float e3[3];
    float c1[3];
    float c2[3];

    e1[0] = v1[0] - v0[0];
    e1[1] = v1[1] - v0[1];
    e1[2] = v1[2] - v0[2];

    e2[0] = v2[0] - v0[0];
    e2[1] = v2[1] - v0[1];
    e2[2] = v2[2] - v0[2];

    e3[0] = v3[0] - v0[0];
    e3[1] = v3[1] - v0[1];
    e3[2] = v3[2] - v0[2];

    c1[0] = e1[1] * e2[2] - e1[2] * e2[1];
    c1[1] = e1[2] * e2[0] - e1[0] * e2[2];
    c1[2] = e1[0] * e2[1] - e1[1] * e2[0];

    c2[0] = e3[1] * e2[2] - e3[2] * e2[1];
    c2[1] = e3[2] * e2[0] - e3[0] * e2[2];
    c2[2] = e3[0] * e2[1] - e3[1] * e2[0];

    float area = 0;
    area += 0.5f * sqrt(c1[0] * c1[0] + c1[1] * c1[1] + c1[2] * c1[2]);
    area += 0.5f * sqrt(c2[0] * c2[0] + c2[1] * c2[1] + c2[2] * c2[2]);

    return area;
}

void Visualizer::wuDrawStreamsurface()
{
    GLfloat positionX;
    GLfloat positionY;
    GLfloat positionZ;
    GLfloat prevPositionX;
    GLfloat prevPositionY;
    GLfloat prevPositionZ;

    GLfloat scale = 1.3f;
    GLfloat h = scale / maxN;

    int x, y, z;
    x = solver->airconditioners[0].pos[0];
    y = solver->airconditioners[0].pos[1];
    z = solver->airconditioners[0].pos[2];

    float seeds[40][3];
    int numSeed = 0;
    for (int i = -2; i <= 2; i++) {
        seeds[numSeed][0] = x - 4;
        seeds[numSeed][1] = y - 1;
        seeds[numSeed][2] = z + i;
        numSeed++;
    }
    for (int i = -2; i <= 2; i++) {
        seeds[numSeed][0] = x + 4;
        seeds[numSeed][1] = y - 1;
        seeds[numSeed][2] = z + i;
        numSeed++;
    }
    for (int i = -2; i <= 2; i++) {
        seeds[numSeed][0] = x + i;
        seeds[numSeed][1] = y - 1;
        seeds[numSeed][2] = z - 4;
        numSeed++;
    }
    for (int i = -2; i <= 2; i++) {
        seeds[numSeed][0] = x + i;
        seeds[numSeed][1] = y - 1;
        seeds[numSeed][2] = z + 4;
        numSeed++;
    }

    if(false)
    {

        for (int i = -2; i <= 2; i++) {
            seeds[numSeed][0] = x - 7;
            seeds[numSeed][1] = y - 2;
            seeds[numSeed][2] = z + i;
            numSeed++;
        }
        for (int i = -2; i <= 2; i++) {
            seeds[numSeed][0] = x + 7;
            seeds[numSeed][1] = y - 2;
            seeds[numSeed][2] = z + i;
            numSeed++;
        }
        for (int i = -2; i <= 2; i++) {
            seeds[numSeed][0] = x + i;
            seeds[numSeed][1] = y - 2;
            seeds[numSeed][2] = z - 7;
            numSeed++;
        }
        for (int i = -2; i <= 2; i++) {
            seeds[numSeed][0] = x + i;
            seeds[numSeed][1] = y - 2;
            seeds[numSeed][2] = z + 7;
            numSeed++;
        }
    }

    for (int i = 0; i < numSeed; i++) {
        seeds[i][0] += 0.5f;
        seeds[i][2] += 0.5f;
    }


    float particlePos[40][stLen][3];

    for (int i = 0; i < numSeed; i++) {

        particlePos[i][0][0] = seeds[i][0];
        particlePos[i][0][1] = seeds[i][1];
        particlePos[i][0][2] = seeds[i][2];

        for (int j = 1; j < stLen; j++) {

            float vx, vy, vz;
            wuGetVelocityRK4(particlePos[i][j - 1][0], particlePos[i][j - 1][1], particlePos[i][j - 1][2], vx, vy, vz);
            particlePos[i][j][0] = particlePos[i][j - 1][0] + particleSpeed * vx;
            particlePos[i][j][1] = particlePos[i][j - 1][1] + particleSpeed * vy;
            particlePos[i][j][2] = particlePos[i][j - 1][2] + particleSpeed * vz;

            wu_no_stick(particlePos[i][j][0], particlePos[i][j][1], particlePos[i][j][2]);

        }

    }


    glColor4f(1.f, 1.f, 1.f, 0.3f);
    glBegin(GL_QUADS);

    for (int batch = 0; batch < 4; batch++) {
        float areaZero = 1;

        for (int j = 1; j < stLen; j++) {
            float area = 1;

            //literal area
            //for (int hole = 0; hole < 4; hole++) {
            //    int i = 5 * batch + hole;
            //    
            //    area += areaQuad(&particlePos[i][j][0], &particlePos[i + 1][j][0], &particlePos[i + 1][j - 1][0], &particlePos[i][j - 1][0]);
            //}

            //distance
            //float d[3];
            //d[0] = particlePos[5 * batch + 4][j][0] - particlePos[5 * batch][j][0];
            //d[1] = particlePos[5 * batch + 4][j][1] - particlePos[5 * batch][j][1];
            //d[2] = particlePos[5 * batch + 4][j][2] - particlePos[5 * batch][j][2];
            //area= sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
            
            area = 1 - (j - 1.0) / (stLen - 2);

            if (j == 1) {
                areaZero = area;
            }

            float t = wuGetDensity(particlePos[5 * batch][j][0], particlePos[5 * batch][j][1], particlePos[5 * batch][j][2]);
            //t = (t - minTem) / (maxTem - minTem);
            t = (t - 15) / (maxTem - 15);

            if (t < 0) t = 0;
            if (t > 1) t = 1;

            //float opacity = areaZero / area - 0.2f;
            // 
            float opacity = area > 0 ? area : 0;
            //float opacity = 2.f - 1.f * area / areaZero;

            if (t < 0.25f) {
                glColor4f(0, 4 * t, 1.f, opacity);
            }
            else if (t < 0.5) {
                glColor4f(0, 1.f, 2.f - 4 * t, opacity);
            }
            else if (t < 0.75) {
                glColor4f(4 * t - 2, 1.f, 0, opacity);
            }
            else {
                glColor4f(1.f, 4 - 4 * t, 0, opacity);
            }
            //ignore coloring
            glColor4f(1.f, 1.f, 1.f, opacity);

            for (int hole = 0; hole < 4; hole++) {
                int i = 5 * batch + hole;

                glVertex3f(particlePos[i][j][0] * h - scale / 2, particlePos[i][j][1] * h, particlePos[i][j][2] * h - scale / 2);
                glVertex3f(particlePos[i + 1][j][0] * h - scale / 2, particlePos[i + 1][j][1] * h, particlePos[i + 1][j][2] * h - scale / 2);
                glVertex3f(particlePos[i + 1][j - 1][0] * h - scale / 2, particlePos[i + 1][j - 1][1] * h, particlePos[i + 1][j - 1][2] * h - scale / 2);
                glVertex3f(particlePos[i][j - 1][0] * h - scale / 2, particlePos[i][j - 1][1] * h, particlePos[i][j - 1][2] * h - scale / 2);
            }
        }
    }

    glEnd();
}

#define maxNodeCnt 70
struct surface {
    GLfloat gridPos[maxNodeCnt][stLen][3];
    int nextIdx[maxNodeCnt][stLen];
    //float areaZero[maxNodeCnt][stLen];
    float areaZero[maxNodeCnt];
    int pred[maxNodeCnt][stLen];
    int succ[maxNodeCnt][stLen];
    bool valid[maxNodeCnt][stLen];
    float opacity[maxNodeCnt][stLen];
};

void Visualizer::wuDrawStreamsurfaceWithSplit()
{
    GLfloat scale = 1.3f;
    GLfloat h = scale / maxN;


    static surface surfaces[4];
    int surfaceNum;
    bool inited = false;
    AirConditioner& ac = solver->airconditioners[0];

    if (!inited) {
        int x, y, z;
        x = ac.pos[0];
        y = ac.pos[1];
        z = ac.pos[2];

        int noteCnt = 0;

        if (ac.acType == ACTYPE::Ceiling) {
            surfaceNum = 4;

            noteCnt = 16;
            int interval = (maxNodeCnt-1) / noteCnt; //=4
            for (int i = 0; i <= noteCnt; i++) {
                int ii = interval * i;
                surfaces[0].gridPos[ii][0][0] = x - 4 + 0.5f;// -0.2f * sin((i + 2 * solver->getSimulationTime()) * 3.141592f / 4.0f);
                surfaces[0].gridPos[ii][0][1] = y - 0.2f;
                surfaces[0].gridPos[ii][0][2] = z + i * 4.0f / noteCnt - 2 + 0.5f;

                surfaces[1].gridPos[ii][0][0] = x + 4 + 0.5f;// +0.2f * sin((i - 2 * solver->getSimulationTime()) * 3.141592f / 4.0f);
                surfaces[1].gridPos[ii][0][1] = y - 0.2f;
                surfaces[1].gridPos[ii][0][2] = z + i * 4.0f / noteCnt - 2 + 0.5f;

                surfaces[2].gridPos[ii][0][0] = x + i * 4.0f / noteCnt - 2 + 0.5f;
                surfaces[2].gridPos[ii][0][1] = y - 0.2f;
                surfaces[2].gridPos[ii][0][2] = z - 4 + 0.5f;// -0.2f * sin((i - 2 * solver->getSimulationTime()) * 3.141592f / 4.0f);

                surfaces[3].gridPos[ii][0][0] = x + i * 4.0f / noteCnt - 2 + 0.5f;
                surfaces[3].gridPos[ii][0][1] = y - 0.2f;
                surfaces[3].gridPos[ii][0][2] = z + 4 + 0.5f;// +0.2f * sin((i + 2 * solver->getSimulationTime()) * 3.141592f / 4.0f);

                surfaces[0].nextIdx[ii][0] = ii + interval;
                surfaces[1].nextIdx[ii][0] = ii + interval;
                surfaces[2].nextIdx[ii][0] = ii + interval;
                surfaces[3].nextIdx[ii][0] = ii + interval;

            }
            surfaces[0].nextIdx[interval * noteCnt][0] = -1;
            surfaces[1].nextIdx[interval * noteCnt][0] = -1;
            surfaces[2].nextIdx[interval * noteCnt][0] = -1;
            surfaces[3].nextIdx[interval * noteCnt][0] = -1;

            for (int s = 0; s < 4; s++) {
                for (int i = 0; i < 16; i++) {

                    surfaces[s].areaZero[interval * i] = Dist(surfaces[s].gridPos[interval * i][0], surfaces[s].gridPos[interval * i + interval][0]);
                    surfaces[s].valid[interval * i][0] = true;
                    surfaces[s].opacity[interval * i][0] = 1.0f -0.2f + 0.2f * sin((i + 2 * solver->getSimulationTime()) * 3.141592f / 4.0f);

                }
                surfaces[s].opacity[interval * 16][0] = 1.0f;
            }
            
        }
        else if (ac.acType == ACTYPE::Tower)
        {
            surfaceNum = 4;
            //printf("AC TYPE: Tower");
            //printf("\n direction:%d", solver->airconditioners[0].direction);
            if (ac.direction % 2 == 0)
            {
                int horX = 0, horZ = 0; //right-end direction
                if (ac.direction == 0) {
                    horX = 0;
                    horZ = 1;
                }
                else if (ac.direction == 2) {
                    horX = 1;
                    horZ = 0;
                }
                else if (ac.direction == 4) {
                    horX = 0;
                    horZ = -1;
                }
                else if (ac.direction == 6) {
                    horX = -1;
                    horZ = 0;
                }

                int interval = 16;
                for (int i = 0; i < 5; i++) {
                    int ii = interval * i;

                    surfaces[0].gridPos[ii][0][0] = x - 2 * horX + 0.5f;
                    surfaces[0].gridPos[ii][0][1] = y + 5 + i + 0.5f;
                    surfaces[0].gridPos[ii][0][2] = z - 2 * horZ + 0.5f;

                    surfaces[1].gridPos[ii][0][0] = x - 2 * horX + 0.5f;
                    surfaces[1].gridPos[ii][0][1] = y + 10 + i + 0.5f;
                    surfaces[1].gridPos[ii][0][2] = z - 2 * horZ + 0.5f;

                    surfaces[2].gridPos[ii][0][0] = x + 1 * horX + 0.5f;
                    surfaces[2].gridPos[ii][0][1] = y + 5 + i + 0.5f;
                    surfaces[2].gridPos[ii][0][2] = z + 1 * horZ + 0.5f;

                    surfaces[3].gridPos[ii][0][0] = x + 1 * horX + 0.5f;
                    surfaces[3].gridPos[ii][0][1] = y + 10 + i + 0.5f;
                    surfaces[3].gridPos[ii][0][2] = z + 1 * horZ + 0.5f;

                    for (int s = 0; s < 4; s++){
                        surfaces[s].nextIdx[ii][0] = (i == 4) ? -1 : ii + interval;
                    }
                }

                for (int s = 0; s < 4; s++) {
                    for (int i = 0; i < 4; i++) {
                        int ii = interval * i;
                        surfaces[s].areaZero[ii] = Dist(surfaces[s].gridPos[ii][0], surfaces[s].gridPos[ii + interval][0]);
                        surfaces[s].valid[ii][0] = true;
                        surfaces[s].opacity[ii][0] = 1.0f;
                    }
                    surfaces[s].opacity[interval * 4][0] = 1.0f;
                }

                if(false)
                {
                    float center[3];
                    center[0] = x - 0.5f * horX + 0.5f;
                    center[1] = y + 17.f;
                    center[2] = z - 0.5f * horZ + 0.5f;
                    float hor[3];
                    hor[0] = -1.5f * horX;
                    hor[1] = 0;
                    hor[2] = -1.5f * horZ;

                    int numRing = 16;
                    for (int i = 0; i <= numRing; i++) {
                        int ii = 4 * i;
                        const float theta = i * 2 * 3.141592f / numRing;
                        surfaces[4].gridPos[ii][0][0] = center[0] + 1.0f * (hor[0] * cos(theta));
                        surfaces[4].gridPos[ii][0][1] = center[1] + 1.0f * 1.5f * sin(theta);
                        surfaces[4].gridPos[ii][0][2] = center[2] + 1.0f * (hor[2] * cos(theta));

                        surfaces[4].nextIdx[ii][0] = ii + 4;
                    }
                    for (int i = 0; i < numRing; i++) {
                        int ii = 4 * i;
                        surfaces[4].areaZero[ii] = Dist(surfaces[4].gridPos[ii][0], surfaces[4].gridPos[ii + 4][0]);
                        surfaces[4].valid[ii][0] = true;
                        surfaces[4].opacity[ii][0] = 1.0f;
                    }
                    surfaces[4].nextIdx[64][0] = -1;
                    surfaces[4].opacity[64][0] = 1.0f;
                }

            }
            else
            {   //diagonal placement

                int horX = 0, horZ = 0; //right-end direction
                if (ac.direction == 1) {
                    horX = 1;
                    horZ = 1;
                }
                else if (ac.direction == 3) {
                    horX = 1;
                    horZ = -1;
                }
                else if (ac.direction == 5) {
                    horX = -1;
                    horZ = -1;
                }
                else if (ac.direction == 7) {
                    horX = -1;
                    horZ = 1;
                }

                int interval = 16;
                for (int i = 0; i < 5; i++) {
                    int ii = interval * i;

                    surfaces[0].gridPos[ii][0][0] = x - horX + 0.5f;
                    surfaces[0].gridPos[ii][0][1] = y + 5 + i + 0.5f;
                    surfaces[0].gridPos[ii][0][2] = z - horZ + 0.5f;

                    surfaces[1].gridPos[ii][0][0] = x - horX + 0.5f;
                    surfaces[1].gridPos[ii][0][1] = y + 10 + i + 0.5f;
                    surfaces[1].gridPos[ii][0][2] = z - horZ + 0.5f;

                    surfaces[2].gridPos[ii][0][0] = x + horX + 0.5f;
                    surfaces[2].gridPos[ii][0][1] = y + 5 + i + 0.5f;
                    surfaces[2].gridPos[ii][0][2] = z + horZ + 0.5f;

                    surfaces[3].gridPos[ii][0][0] = x + horX + 0.5f;
                    surfaces[3].gridPos[ii][0][1] = y + 10 + i + 0.5f;
                    surfaces[3].gridPos[ii][0][2] = z + horZ + 0.5f;

                    for (int s = 0; s < 4; s++) {
                        surfaces[s].nextIdx[ii][0] = (i == 4) ? -1 : ii + interval;
                    }
                }

                for (int s = 0; s < 4; s++) {
                    for (int i = 0; i < 4; i++) {
                        int ii = interval * i;
                        surfaces[s].areaZero[ii] = Dist(surfaces[s].gridPos[ii][0], surfaces[s].gridPos[ii + interval][0]);
                        surfaces[s].valid[ii][0] = true;
                        surfaces[s].opacity[ii][0] = 1.0f;
                    }
                    surfaces[s].opacity[interval * 4][0] = 1.0f;
                }

                if(false)
                {
                    float center[3];
                    center[0] = x + 0.5f;
                    center[1] = y + 17.f;
                    center[2] = z + 0.5f;
                    float hor[3];
                    hor[0] = -1.5f / sqrt(2) * horX;
                    hor[1] = 0;
                    hor[2] = -1.5f / sqrt(2) * horZ;

                    int numRing = 16;
                    for (int i = 0; i <= numRing; i++) {
                        int ii = 4 * i;
                        const float theta = i * 2 * 3.141592f / numRing;
                        surfaces[4].gridPos[ii][0][0] = center[0] + 1.0f * (hor[0] * cos(theta));
                        surfaces[4].gridPos[ii][0][1] = center[1] + 1.0f * 1.5f * sin(theta);
                        surfaces[4].gridPos[ii][0][2] = center[2] + 1.0f * (hor[2] * cos(theta));

                        surfaces[4].nextIdx[ii][0] = ii + 4;
                    }
                    for (int i = 0; i < numRing; i++) {
                        int ii = 4 * i;
                        surfaces[4].areaZero[ii] = Dist(surfaces[4].gridPos[ii][0], surfaces[4].gridPos[ii + 4][0]);
                        surfaces[4].valid[ii][0] = true;
                        surfaces[4].opacity[ii][0] = 1.0f;
                    }
                    surfaces[4].nextIdx[64][0] = -1;
                    surfaces[4].opacity[64][0] = 1.0f;
                }

            }
        }
        else if (ac.acType == ACTYPE::Wall)
        {
            surfaceNum = 1;


            int interval = 9;
            for (int s = 0; s < surfaceNum; s++) {
                for (int i = 0; i <= 7; i++) { //7=2*widthSide+1
                    int ii = interval * i;
                    if (ac.direction % 4 == 0) {
                        surfaces[s].gridPos[ii][0][0] = x + 0.5f;
                        surfaces[s].gridPos[ii][0][1] = y + 3*(s -1 ) + 0.5f;
                        surfaces[s].gridPos[ii][0][2] = z + (i - 3);
                    }
                    else {
                        surfaces[s].gridPos[ii][0][0] = x + (i - 3);
                        surfaces[s].gridPos[ii][0][1] = y + 3*(s-1) + 0.5f;
                        surfaces[s].gridPos[ii][0][2] = z + 0.5f;
                    }
                    surfaces[s].nextIdx[ii][0] = ii + interval;
                }
                surfaces[s].nextIdx[interval * 6][0] = -1;
            }


            for (int s = 0; s < 3; s++) {
                for (int i = 0; i < 7; i++) {

                    surfaces[s].areaZero[interval * i] = Dist(surfaces[s].gridPos[interval * i][0], surfaces[s].gridPos[interval * i + interval][0]);
                    surfaces[s].valid[interval * i][0] = true;
                    surfaces[s].opacity[interval * i][0] = 1.0f;

                }
                surfaces[s].opacity[interval * 7][0] = 1.0f;
            }


        }
        else if (ac.acType == ACTYPE::Stand)
        {
            surfaceNum = 1;
            int nodeCnt = 10;
            int interval = 69 / nodeCnt;
            for (int s = 0; s < surfaceNum; s++) {
                for (int i = 0; i <= nodeCnt; i++) {
                    int ii = interval * i;
                    if (ac.direction % 4 == 0) {
                        surfaces[s].gridPos[ii][0][0] = x + 0.5f;// -0.2f * sin((i + 2 * solver->getSimulationTime()) * 3.141592f / 4.0f);
                        surfaces[s].gridPos[ii][0][1] = y + 1*(s) + 0.5f;
                        surfaces[s].gridPos[ii][0][2] = z + (0.5f * i - 2);
                        surfaces[s].nextIdx[ii][0] = ii + interval;

                    }
                    else {
                        surfaces[s].gridPos[ii][0][0] = x + (0.5f * i - 2);
                        surfaces[s].gridPos[ii][0][1] = y + 1*(s) + 0.5f;
                        surfaces[s].gridPos[ii][0][2] = z + 0.5f;
                        surfaces[s].nextIdx[ii][0] = ii + interval;


                    }
                }
                surfaces[s].nextIdx[interval * nodeCnt][0] = -1;
            }

            for (int s = 0; s < surfaceNum; s++) {
                for (int i = 0; i < nodeCnt; i++) {

                    surfaces[s].areaZero[interval * i] = Dist(surfaces[s].gridPos[interval * i][0], surfaces[s].gridPos[interval * i + interval][0]);
                    surfaces[s].valid[interval * i][0] = true;
                    surfaces[s].opacity[interval * i][0] = 1.0f;

                    surfaces[s].areaZero[interval * i] *= 1 - 0.1f + 0.1f * sin((2*i + 2 * solver->getSimulationTime()) * 3.141592f / 4.0f);
                }
                surfaces[s].opacity[interval * nodeCnt][0] = 1.0f;
            }

        }

        inited = true;

    }

    const float minWidth = 0.04f;
    const float maxWidth = 0.4f;
    float alpha = 2.0f;
    float beta = 1.3f;
    float delta = 0.8f;
    float zeta = 0.8f;
    float gamma = 0.02f;
    float K = 0.002f;
    float s_shape = 0.7f;
    //int Len = stLen;

    glColor4f(1.f, 1.f, 1.f, 0.3f);
    glBegin(GL_TRIANGLES);
    for (int s = 0; s < surfaceNum; s++) {
        for (int j = 1; j < stLen; j++) {


            //advection
            int i_1 = 0;
            int i = 0;
            for (i = 0; i != -1; i_1 = i, i = surfaces[s].nextIdx[i][j - 1]) {

                surfaces[s].pred[i][j] = i;
                surfaces[s].succ[i][j - 1] = i;
                surfaces[s].valid[i][j] = surfaces[s].valid[i][j - 1];


                if (surfaces[s].valid[i][j - 1] || surfaces[s].valid[i_1][j - 1])
                {
                    float vx, vy, vz;
                    wuGetVelocityRK4(surfaces[s].gridPos[i][j - 1][0], surfaces[s].gridPos[i][j - 1][1], surfaces[s].gridPos[i][j - 1][2], vx, vy, vz);
                    
                    {   //normalize
                        float mag = sqrt(vx * vx + vy * vy + vz * vz);
                        if (mag <= 1e-5) {
                            surfaces[s].valid[i][j] = false;
                            vx = vy = vz = 0;
                        }
                        else {
                            vx /= mag * 4;
                            vy /= mag * 4;
                            vz /= mag * 4;
                        }
                    }
                    surfaces[s].gridPos[i][j][0] = surfaces[s].gridPos[i][j - 1][0] + 0.3f * particleSpeed * vx;
                    surfaces[s].gridPos[i][j][1] = surfaces[s].gridPos[i][j - 1][1] + 0.3f * particleSpeed * vy;
                    surfaces[s].gridPos[i][j][2] = surfaces[s].gridPos[i][j - 1][2] + 0.3f * particleSpeed * vz;

                    wu_no_stick(surfaces[s].gridPos[i][j][0], surfaces[s].gridPos[i][j][1], surfaces[s].gridPos[i][j][2]);
                }

                if (i != 0)
                    surfaces[s].nextIdx[i_1][j] = i;
            }
            surfaces[s].nextIdx[i_1][j] = -1;

            if(false)
            {   //splitting
                i_1 = i = 0;
                for (int i1 = surfaces[s].nextIdx[0][j]; i1 != -1; i_1 = i, i = i1, i1 = surfaces[s].nextIdx[i1][j]) {

                    if (i1 - i <= 1)
                        continue;

                    if (!surfaces[s].valid[i][j])
                        continue;

                    int i_new = (i + i1) / 2;
                    int i2 = surfaces[s].nextIdx[i1][j];

                    bool splitFlag = false;

                    //alpha
                    //if (Dist(surfaces[s].gridPos[i][j], surfaces[s].gridPos[i1][j]) > alpha * surfaces[s].areaZero[i][j]) {
                    if (Dist(surfaces[s].gridPos[i][j], surfaces[s].gridPos[i1][j]) > maxWidth) {
                        splitFlag = true;
                        //printf("\n split at (%d,%d,%d,%d) / criteria 1", s, j, i, i1);
                    }

                    ////beta
                    if (i2 != -1 && surfaces[s].valid[i_1][j] && surfaces[s].valid[i][j] && surfaces[s].valid[i1][j] &&
                        Curv(surfaces[s].gridPos[i_1][j], surfaces[s].gridPos[i][j], surfaces[s].gridPos[i1][j])
                        + Curv(surfaces[s].gridPos[i][j], surfaces[s].gridPos[i1][j], surfaces[s].gridPos[i2][j]) > beta)
                    {
                        splitFlag = true;
                        //printf("\n split at (%d,%d,%d,%d) / criteria 2", s, j, i, i1);
                    }

                    //if (j % 100 == 0)
                    //    splitFlag = true;
                    //else
                    //    splitFlag = false;
                    if (!splitFlag)
                        continue;


                    //if(true){
                    if (i == 0 || i2 == -1) {
                        surfaces[s].gridPos[i_new][j][0] = (surfaces[s].gridPos[i][j][0] + surfaces[s].gridPos[i1][j][0]) / 2.0f;
                        surfaces[s].gridPos[i_new][j][1] = (surfaces[s].gridPos[i][j][1] + surfaces[s].gridPos[i1][j][1]) / 2.0f;
                        surfaces[s].gridPos[i_new][j][2] = (surfaces[s].gridPos[i][j][2] + surfaces[s].gridPos[i1][j][2]) / 2.0f;
                        //wu_no_stick(surfaces[s].gridPos[i_new][j][0], surfaces[s].gridPos[i_new][j][1], surfaces[s].gridPos[i_new][j][2]);
                    }
                    else {
                        surfaces[s].gridPos[i_new][j][0] = (surfaces[s].gridPos[i][j][0] + surfaces[s].gridPos[i1][j][0]) * 9.0f / 16.0f - (surfaces[s].gridPos[i_1][j][0] + surfaces[s].gridPos[i2][j][0]) / 16.0f;
                        surfaces[s].gridPos[i_new][j][1] = (surfaces[s].gridPos[i][j][1] + surfaces[s].gridPos[i1][j][1]) * 9.0f / 16.0f - (surfaces[s].gridPos[i_1][j][1] + surfaces[s].gridPos[i2][j][1]) / 16.0f;
                        surfaces[s].gridPos[i_new][j][2] = (surfaces[s].gridPos[i][j][2] + surfaces[s].gridPos[i1][j][2]) * 9.0f / 16.0f - (surfaces[s].gridPos[i_1][j][2] + surfaces[s].gridPos[i2][j][2]) / 16.0f;
                        //wu_no_stick(surfaces[s].gridPos[i_new][j][0], surfaces[s].gridPos[i_new][j][1], surfaces[s].gridPos[i_new][j][2]);
                    }
                    float d1 = Dist(surfaces[s].gridPos[i][j], surfaces[s].gridPos[i_new][j]);
                    float d2 = Dist(surfaces[s].gridPos[i_new][j], surfaces[s].gridPos[i1][j]);

                    //surfaces[s].areaZero[i_new][j] = surfaces[s].areaZero[i][j] * d1 / (d1 + d2); // 0.5f
                    //surfaces[s].areaZero[i][j] = surfaces[s].areaZero[i][j] - surfaces[s].areaZero[i_new][j];
                    surfaces[s].areaZero[i_new] = surfaces[s].areaZero[i] * d1 / (d1 + d2); // 0.5f
                    surfaces[s].areaZero[i] = surfaces[s].areaZero[i] - surfaces[s].areaZero[i_new];

                    //printf("\n at (%d,%d,%d,%d): splitt %f into %f + %f", s, i, i_new, j, surfaces[s].areaZero[i]+ surfaces[s].areaZero[i_new], surfaces[s].areaZero[i], surfaces[s].areaZero[i_new]);

                    surfaces[s].nextIdx[i][j] = i_new;
                    surfaces[s].pred[i_new][j] = i1;
                    surfaces[s].valid[i_new][j] = true;
                    surfaces[s].nextIdx[i_new][j] = i1;
                }
            }

            if(false)
            {   //removing
                i_1 = 0;
                i = surfaces[s].nextIdx[0][j];
                if (i == -1) continue;
                for (int i1 = surfaces[s].nextIdx[i][j]; i1 != -1; i = i1, i1 = surfaces[s].nextIdx[i1][j]) {

                    int i2 = surfaces[s].nextIdx[i1][j];
                    if (i2 == -1)
                        break;

                    if (!surfaces[s].valid[i_1][j] || !surfaces[s].valid[i][j] || !surfaces[s].valid[i1][j]) {
                        i_1 = i;
                        continue;
                    }

                    bool mergeFlag = false;
                    //float sumArea = surfaces[s].areaZero[i_1][j] + surfaces[s].areaZero[i][j];
                    float sumArea = surfaces[s].areaZero[i_1] + surfaces[s].areaZero[i];

                    //delta and zeta
                    //if (Dist(surfaces[s].gridPos[i_1][j], surfaces[s].gridPos[i][j]) + Dist(surfaces[s].gridPos[i1][j], surfaces[s].gridPos[i][j]) < delta * sumArea)
                    if (Dist(surfaces[s].gridPos[i_1][j], surfaces[s].gridPos[i][j]) + Dist(surfaces[s].gridPos[i1][j], surfaces[s].gridPos[i][j]) < 2 * minWidth)
                        //if (Curv(surfaces[s].gridPos[i_1][j], surfaces[s].gridPos[i][j], surfaces[s].gridPos[i1][j]) + Curv(surfaces[s].gridPos[i][j], surfaces[s].gridPos[i1][j], surfaces[s].gridPos[i2][j]) < zeta)
                    {
                        mergeFlag = true;
                        //printf("\n merge at (%d,%d,%d,%d)", s, j, i, i1);
                    }

                    //mergeFlag = false;
                    if (!mergeFlag) {
                        i_1 = i;
                        continue;
                    }
                    //printf("\n at (%d,%d,%d,%d): merged %f + %f into %f", s, i_1, i, j, surfaces[s].areaZero[i_1][j], surfaces[s].areaZero[i][j], sumArea);
                    int pre = surfaces[s].pred[i][j];
                    while (surfaces[s].succ[pre][j - 1] == i) {
                        surfaces[s].succ[pre][j - 1] = i_1;
                        pre = surfaces[s].nextIdx[pre][j - 1];
                    }
                    //surfaces[s].areaZero[i_1][j] = sumArea;
                    //surfaces[s].areaZero[i][j] = 0;
                    surfaces[s].areaZero[i_1] = sumArea;
                    surfaces[s].areaZero[i] = 0;

                    surfaces[s].valid[i][j] = false;
                    surfaces[s].nextIdx[i_1][j] = i1;

                    i_1 = i;

                }
            }

            {   //validity check
                i = 0;
                for (int i1 = surfaces[s].nextIdx[0][j]; i1 != -1; i = i1, i1 = surfaces[s].nextIdx[i1][j]) {
                    int _i = surfaces[s].pred[i][j];
                    int _i1 = surfaces[s].pred[i1][j];

                    surfaces[s].opacity[i][j] = 0;

                    if (_i == _i1)
                        continue;
                    if (!surfaces[s].valid[i][j])
                        continue;

                    float d1 = Dist(surfaces[s].gridPos[i][j], surfaces[s].gridPos[i1][j]);
                    //float d2 = Dist(surfaces[s].gridPos[_i][j - 1], surfaces[s].gridPos[_i1][j - 1]);
                    //float d3 = Dist(surfaces[s].gridPos[i][j], surfaces[s].gridPos[_i][j - 1]);

                    ////gamma
                    //if (d1 - d2 > gamma * d3) {
                    //    surfaces[s].valid[i][j] = false;
                    //    //printf("\n invalidation at (%d,%d,%d,%d)", s, j, i, i1);
                    //}


                    //if (d1 > 3 * surfaces[s].areaZero[i])
                    //    surfaces[s].valid[i][j] = false;

                    if(surfaces[s].gridPos[i][j][1] > solver->dim.y-2.7f)
                    {
                        surfaces[s].valid[i][j] = false;
                        //printf("\n invalidation at (%d,%d,%d,%d)", s, j, i, i1);
                    }

                    //if(false)
                    {
                        float D[3];
                        float V[3];
                        for (int q = 0; q < 3; q++) {
                            D[q] = surfaces[s].gridPos[_i][j - 1][q] - ac.pos[q];
                            V[q] = surfaces[s].gridPos[i][j][q] - surfaces[s].gridPos[_i][j - 1][q];
                        }
                        float DD = sqrt(D[0] * D[0] + D[1] * D[1] + D[2] * D[2]);
                        if (DD < 15) {
                            float VV = sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
                            float DV = D[0] * V[0] + D[1] * V[1] + D[2] * V[2];
                            if (DV / (DD * VV) < 0) {
                                surfaces[s].valid[i][j] = false;
                            }
                        }
                    }

                    //float v_diff[3];
                    //v_diff[0] = (surfaces[s].gridPos[i][j][0] - surfaces[s].gridPos[_i][j - 1][0]) - (surfaces[s].gridPos[i1][j][0] - surfaces[s].gridPos[_i1][j - 1][0]);
                    //v_diff[1] = (surfaces[s].gridPos[i][j][1] - surfaces[s].gridPos[_i][j - 1][1]) - (surfaces[s].gridPos[i1][j][1] - surfaces[s].gridPos[_i1][j - 1][1]);
                    //v_diff[2] = (surfaces[s].gridPos[i][j][2] - surfaces[s].gridPos[_i][j - 1][2]) - (surfaces[s].gridPos[i1][j][2] - surfaces[s].gridPos[_i1][j - 1][2]);

                    //float d2 = v_diff[0] * v_diff[0] + v_diff[1] * v_diff[1] + v_diff[2] * v_diff[2];
                    //if (d2 > 0.05)
                    //{
                    //    surfaces[s].valid[i][j] = false;
                    //    //printf("\n invalidation at (%d,%d,%d,%d)", s, j, i, i1);
                    //}
                }

                surfaces[s].opacity[i][j] = 0;
            }

            {   //opacity calculation
                i = 0;
                i_1 = 0;
                for (int i1 = surfaces[s].nextIdx[0][j]; i1 != -1; i_1 = i, i = i1, i1 = surfaces[s].nextIdx[i1][j]) {

                    //if (!surfaces[s].valid[i][j]) {
                    //    continue;
                    //}

                    float area = 0;
                    //float areaZero = surfaces[s].areaZero[i][j];
                    float areaZero = surfaces[s].areaZero[i];

                    //linear area
                    {
                        area = Dist(surfaces[s].gridPos[i1][j], surfaces[s].gridPos[i][j]);
                    }

                    //quad area
                    if (false) {
                        int _i = surfaces[s].pred[i][j];
                        int _i1 = surfaces[s].pred[i1][j];
                        area = areaQuad(surfaces[s].gridPos[i][j], surfaces[s].gridPos[i1][j], surfaces[s].gridPos[_i1][j - 1], surfaces[s].gridPos[_i][j - 1]);
                        if (area < 0)
                            area = -area;
                    }
                    //area = areaZero;
                    ////////
                    //opacity
                    float opacity = 0;

                    if (true)
                    {
                        if (area <= 1e-9) {
                            opacity = 1.0f;
                        }
                        else {
                            opacity = areaZero / area;
                            //opacity *= 0.8f;
                            if (opacity > 1)
                                opacity = 1;
                            //opacity = opacity * opacity;
                        }
                        opacity *= ((stLen - 1.0f*j) / stLen);
                        
                        if (i < 10)
                            opacity *= i / 10.f;
                        if (i > 64 - 10)
                            opacity *= (64 - i) / 10.f;

                        {   //animation

                            static const float animationSpeed = 40;
                            //static const float phases[16] = { 0.1,6.2,4.3,9.4,4.5,5.6,2.7,6.8,0.9,2.0,8.1,4.2,8.3,1.4,5.5,7.6 };
                            //float p = 2.5f* s + (solver->getFramestep() % 500) / 50.0f - phases[i / 4] - (j % 10000) / 1000.0f;
                            float p = 25 * s + (((int)(solver->getFramestep() * 0.005f * animationSpeed) / 5.0f)) - HashBig(i / 2, 100) - (j % 10000) / 100.0f;
                            while (p < 0)p += 100;
                            while (p >= 100)p -= 100;

                            opacity *= p < 10 ? 1 + (10 - p) * 0.05f : 1;
                        }

                        if(false)
                        {   //red coloring
                            int _i = surfaces[s].pred[i][j];
                            int _i1 = surfaces[s].pred[i1][j];

                            float v_diff[3];
                            v_diff[0] = (surfaces[s].gridPos[i][j][0] - surfaces[s].gridPos[_i][j - 1][0]) - (surfaces[s].gridPos[i1][j][0] - surfaces[s].gridPos[_i1][j - 1][0]);
                            v_diff[1] = (surfaces[s].gridPos[i][j][1] - surfaces[s].gridPos[_i][j - 1][1]) - (surfaces[s].gridPos[i1][j][1] - surfaces[s].gridPos[_i1][j - 1][1]);
                            v_diff[2] = (surfaces[s].gridPos[i][j][2] - surfaces[s].gridPos[_i][j - 1][2]) - (surfaces[s].gridPos[i1][j][2] - surfaces[s].gridPos[_i1][j - 1][2]);

                            float d2 = v_diff[0] * v_diff[0] + v_diff[1] * v_diff[1] + v_diff[2] * v_diff[2];
                            float d3 = Dist(surfaces[s].gridPos[i][j], surfaces[s].gridPos[_i][j - 1]);
                            //if (d2 > d3*gamma)
                            //if (surfaces[s].gridPos[i][j][1] > solver->dim.y - 2.0f)
                            //{
                            //    opacity = 2;
                            //}
                            //if (surfaces[s].opacity[_i][j - 1] == 2)
                            //    opacity = 2;
                        }
                    }
                    else {  //smoke
                        int _i = surfaces[s].pred[i][j];

                        float a_density;
                        float a_curv;

                        float R[3];
                        float V[3];
                        float W[3];
                        float C[3];
                        float cos_gamma;
                        float areaTri;

                        //density

                        R[0] = surfaces[s].gridPos[i][j][0] - translationX;
                        R[1] = surfaces[s].gridPos[i][j][1] - translationY;
                        R[2] = surfaces[s].gridPos[i][j][2] - translationZ;

                        for (int q = 0; q < 3; q++) {
                            V[q] = surfaces[s].gridPos[i1][j][q] - surfaces[s].gridPos[i][j][q];
                            W[q] = surfaces[s].gridPos[_i][j - 1][q] - surfaces[s].gridPos[i][j][q];
                        }

                        C[0] = V[1] * W[2] - V[2] * W[1];
                        C[1] = V[2] * W[0] - V[0] * W[2];
                        C[2] = V[0] * W[1] - V[1] * W[0];

                        areaTri = 0.5f * sqrt(C[0] * C[0] + C[1] * C[1] + C[2] * C[2]);

                        cos_gamma = (C[0] * R[0] + C[1] * R[1] + C[2] * R[2]) / sqrt((C[0] * C[0] + C[1] * C[1] + C[2] * C[2]) * (R[0] * R[0] + R[1] * R[1] + R[2] * R[2]));
                        if (cos_gamma < 0) cos_gamma = -cos_gamma;

                        a_density = K / (cos_gamma * areaTri);
                        if (a_density < 0) a_density = 0;
                        if (a_density > 1) a_density = 1;

                        opacity = a_density;

                        //curv

                        if (i != 0)
                            a_curv = sqrt(Curv(surfaces[s].gridPos[i_1][j], surfaces[s].gridPos[i][j], surfaces[s].gridPos[i1][j]));
                        else
                            a_curv = 0;

                        a_curv = 1 - 2 * a_curv;
                        if (a_curv < 0) a_curv = 0;
                        if (a_curv > 1) a_curv = 1;

                        opacity *= a_curv;
                        //fade

                        opacity *= 1 - (float)j / stLen;
                    }

                    if (surfaces[s].valid[i, j]) {
                        if (surfaces[s].opacity[i][j] == 0)
                            surfaces[s].opacity[i][j] = opacity;
                        else
                            surfaces[s].opacity[i][j] = (surfaces[s].opacity[i][j] + opacity) / 2;
                        surfaces[s].opacity[i1][j] = opacity;
                    }

                }
            }

        }

        //surface rendering
        for (int j = 0; j < stLen; j++) {

            int i = 0;
            for (int i1 = surfaces[s].nextIdx[0][j]; i1 != -1; i = i1, i1 = surfaces[s].nextIdx[i1][j]) {

                int _i1 = surfaces[s].pred[i1][j];
                int i_ = surfaces[s].succ[i][j];

                if (!surfaces[s].valid[i][j]) {
                    //printf("\ninvalid: [%d][%d][%d]", s, i, j);
                    continue;
                }

                //glColor4f(1.f, 1.f, 1.f, 1.f);

                //if(j!=0){
                //    glColor4f(1.f, 1.f, 1.f, surfaces[s].opacity[i][j]);
                //    glVertex3f(surfaces[s].gridPos[i][j][0] * h - scale / 2, surfaces[s].gridPos[i][j][1] * h, surfaces[s].gridPos[i][j][2] * h - scale / 2);
                //    glColor4f(1.f, 1.f, 1.f, surfaces[s].opacity[i1][j]); 
                //    glVertex3f(surfaces[s].gridPos[i1][j][0] * h - scale / 2, surfaces[s].gridPos[i1][j][1] * h, surfaces[s].gridPos[i1][j][2] * h - scale / 2);
                //    glColor4f(1.f, 1.f, 1.f, surfaces[s].opacity[_i1][j - 1]);
                //    glVertex3f(surfaces[s].gridPos[_i1][j - 1][0] * h - scale / 2, surfaces[s].gridPos[_i1][j - 1][1] * h, surfaces[s].gridPos[_i1][j - 1][2] * h - scale / 2);
                //}
                //if (j != stLen - 1) {
                //    glColor4f(1.f, 1.f, 1.f, surfaces[s].opacity[i][j]);
                //    glVertex3f(surfaces[s].gridPos[i][j][0] * h - scale / 2, surfaces[s].gridPos[i][j][1] * h, surfaces[s].gridPos[i][j][2] * h - scale / 2);
                //    glColor4f(1.f, 1.f, 1.f, surfaces[s].opacity[i1][j]);
                //    glVertex3f(surfaces[s].gridPos[i1][j][0] * h - scale / 2, surfaces[s].gridPos[i1][j][1] * h, surfaces[s].gridPos[i1][j][2] * h - scale / 2);
                //    glColor4f(1.f, 1.f, 1.f, surfaces[s].opacity[i_][j + 1]);
                //    glVertex3f(surfaces[s].gridPos[i_][j + 1][0] * h - scale / 2, surfaces[s].gridPos[i_][j + 1][1] * h, surfaces[s].gridPos[i_][j + 1][2] * h - scale / 2);
                //}

                float colori = wuGetDensity(surfaces[s].gridPos[i][j][0], surfaces[s].gridPos[i][j][1], surfaces[s].gridPos[i][j][2]);
                float colori1 = wuGetDensity(surfaces[s].gridPos[i1][j][0], surfaces[s].gridPos[i1][j][1], surfaces[s].gridPos[i1][j][2]);

                if (j != 0) {
                    float color_i1 = wuGetDensity(surfaces[s].gridPos[_i1][j - 1][0], surfaces[s].gridPos[_i1][j - 1][1], surfaces[s].gridPos[_i1][j - 1][2]);
                    wuTemptoGLColor(colori, surfaces[s].opacity[i][j]);
                    glVertex3f(surfaces[s].gridPos[i][j][0] * h - h - scale / 2, surfaces[s].gridPos[i][j][1] * h - h, surfaces[s].gridPos[i][j][2] * h - h - scale / 2);
                    wuTemptoGLColor(colori1, surfaces[s].opacity[i1][j]);
                    glVertex3f(surfaces[s].gridPos[i1][j][0] * h - h - scale / 2, surfaces[s].gridPos[i1][j][1] * h - h, surfaces[s].gridPos[i1][j][2] * h - h - scale / 2);
                    wuTemptoGLColor(color_i1, surfaces[s].opacity[_i1][j - 1]);
                    glVertex3f(surfaces[s].gridPos[_i1][j - 1][0] * h - h - scale / 2, surfaces[s].gridPos[_i1][j - 1][1] * h - h, surfaces[s].gridPos[_i1][j - 1][2] * h - h - scale / 2);
                }
                if (j != stLen - 1) {
                    float colori_ = wuGetDensity(surfaces[s].gridPos[i_][j + 1][0], surfaces[s].gridPos[i_][j + 1][1], surfaces[s].gridPos[i_][j + 1][2]);
                    wuTemptoGLColor(colori, surfaces[s].opacity[i][j]);
                    glVertex3f(surfaces[s].gridPos[i][j][0] * h - h - scale / 2, surfaces[s].gridPos[i][j][1] * h - h, surfaces[s].gridPos[i][j][2] * h - h - scale / 2);
                    wuTemptoGLColor(colori1, surfaces[s].opacity[i1][j]);
                    glVertex3f(surfaces[s].gridPos[i1][j][0] * h - h - scale / 2, surfaces[s].gridPos[i1][j][1] * h - h, surfaces[s].gridPos[i1][j][2] * h - h - scale / 2);
                    wuTemptoGLColor(colori_, surfaces[s].opacity[i_][j + 1]);
                    glVertex3f(surfaces[s].gridPos[i_][j + 1][0] * h - h - scale / 2, surfaces[s].gridPos[i_][j + 1][1] * h - h, surfaces[s].gridPos[i_][j + 1][2] * h - h - scale / 2);
                }
            }
        }
    }


    glEnd();
}

void Visualizer::wuDrawStreakline()
{
    GLfloat positionX;
    GLfloat positionY;
    GLfloat positionZ;
    GLfloat prevPositionX;
    GLfloat prevPositionY;
    GLfloat prevPositionZ;
    GLfloat scale = 1.3f;
    GLfloat h = scale / maxN;


    for (int i = 0; i < visible; i++) {
        if (particles[i].traceLength > 0 && particles[i].traceLength < lineLength) {
            particles[i].particlePos[particles[i].traceLength][0] = particles[i].originalPos[0];
            particles[i].particlePos[particles[i].traceLength][1] = particles[i].originalPos[1];
            particles[i].particlePos[particles[i].traceLength][2] = particles[i].originalPos[2];

            particles[i].traceLength++;
        }
    };

    //int tnow = floor(solver->getSimulationTime());
    //if (tnow > numParticles) tnow = numParticles;
    static int frame = 0;
    static int hole = -2;
    frame++;

    AirConditioner& ac = solver->airconditioners[0];
    int x, y, z;
    x = ac.pos[0];
    y = ac.pos[1];
    z = ac.pos[2];

    if (ac.acType == ACTYPE::Ceiling) {
        if (visible < numParticles - 20 && frame >= 18) {


            hole++;
            if (hole > 2) hole -= 5;

            particles[visible].particlePos[0][0] = particles[visible].originalPos[0] = x - 4 + 0.5f;
            particles[visible].particlePos[0][1] = particles[visible].originalPos[1] = y - 1 + 0.5f;
            particles[visible].particlePos[0][2] = particles[visible].originalPos[2] = z + hole + 0.5f;
            particles[visible].traceLength = 1;
            visible++;

            particles[visible].particlePos[0][0] = particles[visible].originalPos[0] = x + 4 + 0.5f;
            particles[visible].particlePos[0][1] = particles[visible].originalPos[1] = y - 1 + 0.5f;
            particles[visible].particlePos[0][2] = particles[visible].originalPos[2] = z + hole + 0.5f;
            particles[visible].traceLength = 1;
            visible++;

            particles[visible].particlePos[0][0] = particles[visible].originalPos[0] = x + hole + 0.5f;
            particles[visible].particlePos[0][1] = particles[visible].originalPos[1] = y - 1 + 0.5f;
            particles[visible].particlePos[0][2] = particles[visible].originalPos[2] = z - 4 + 0.5f;
            particles[visible].traceLength = 1;
            visible++;

            particles[visible].particlePos[0][0] = particles[visible].originalPos[0] = x + hole + 0.5f;
            particles[visible].particlePos[0][1] = particles[visible].originalPos[1] = y - 1 + 0.5f;
            particles[visible].particlePos[0][2] = particles[visible].originalPos[2] = z + 4 + 0.5f;
            particles[visible].traceLength = 1;
            visible++;

            frame = 0;
        }
    }
    else if (ac.acType == ACTYPE::Tower) {
        if (visible < numParticles - 32 && frame >= 100) {
            if (ac.direction % 2 == 0 && ac.direction < 8) {
                bool isXdir = ac.direction % 4 == 0;
                int flipped = isXdir ? (ac.direction == 0 ? 1 : -1) : (ac.direction == 6 ? 1 : -1);

                //circle
                {
                    for (int i = 0; i < 12; i++) {
                        float angle = i / 12.f * 3.141592f / 180.f;
                        
                        particles[visible].particlePos[0][0] = particles[visible].originalPos[0] = x + 0.5f;
                        particles[visible].particlePos[0][1] = particles[visible].originalPos[1] = 17.5f + 1.5f * cos(i) + 0.5f;
                        particles[visible].particlePos[0][2] = particles[visible].originalPos[2] = z - flipped * 0.5f + 1.5f * sin(i) + 0.5f;
                        particles[visible].traceLength = 1;
                        visible++;
                    }
                }

                //guard
                for (int i = 6; i <= 15; i++) {


                    particles[visible].particlePos[0][0] = particles[visible].originalPos[0] = x + 0.5f;
                    particles[visible].particlePos[0][1] = particles[visible].originalPos[1] = i + 0.5f;
                    particles[visible].particlePos[0][2] = particles[visible].originalPos[2] = z - 2 * flipped + 0.5f;
                    particles[visible].traceLength = 1;
                    visible++;

                    particles[visible].particlePos[0][0] = particles[visible].originalPos[0] = x + 0.5f;
                    particles[visible].particlePos[0][1] = particles[visible].originalPos[1] = i + 0.5f;
                    particles[visible].particlePos[0][2] = particles[visible].originalPos[2] = z + 1 * flipped + 0.5f;
                    particles[visible].traceLength = 1;
                    visible++;

                    //particles[visible].particlePos[0][0] = particles[visible].originalPos[0] = x + hole + 0.5f;
                    //particles[visible].particlePos[0][1] = particles[visible].originalPos[1] = i + 0.5f;
                    //particles[visible].particlePos[0][2] = particles[visible].originalPos[2] = z - 0.7f * flipped + 0.5f;
                    //particles[visible].traceLength = 1;
                    //visible++;

                    //particles[visible].particlePos[0][0] = particles[visible].originalPos[0] = x + hole + 0.5f;
                    //particles[visible].particlePos[0][1] = particles[visible].originalPos[1] = i + 0.5f;
                    //particles[visible].particlePos[0][2] = particles[visible].originalPos[2] = z + 0 * flipped + 0.5f;
                    //particles[visible].traceLength = 1;
                    //visible++;
                }
            }

            //reset
            frame = 0;

        }

    }

    glColor3f(1.f, 1.f, 1.f);
    glBegin(GL_LINES);
    //glBegin(GL_POINTS);

    for (int i = 0; i < visible; i++) {
        if (particles[i].age > lifespan) continue;

        for (int j = 1; j < particles[i].traceLength; j++) {

            positionX = particles[i].particlePos[j][0] * h - h - scale / 2.f;
            positionY = particles[i].particlePos[j][1] * h - h;
            positionZ = particles[i].particlePos[j][2] * h - h - scale / 2.f;

            prevPositionX = particles[i].particlePos[j - 1][0] * h - h - scale / 2.f;
            prevPositionY = particles[i].particlePos[j - 1][1] * h - h;
            prevPositionZ = particles[i].particlePos[j - 1][2] * h - h - scale / 2.f;

            glVertex3f(positionX, positionY, positionZ);
            glVertex3f(prevPositionX, prevPositionY, prevPositionZ);

            //glVertex3f(positionX, positionY, positionZ);
        }
    }
    glEnd();

}

void Visualizer::wuDrawStreamsurfaceWithLines()
{
    wuDrawStreamsurfaceWithSplit();
    wuDrawStreamline(true);
}

void Visualizer::wuDrawStreaksurface()
{

    GLfloat positionX;
    GLfloat positionY;
    GLfloat positionZ;
    GLfloat prevPositionX;
    GLfloat prevPositionY;
    GLfloat prevPositionZ;
    GLfloat scale = 1.3f;
    GLfloat h = scale / maxN;


    for (int i = 0; i < visible; i++) {
        if (particles[i].traceLength > 0 && particles[i].traceLength < lineLength) {
            particles[i].particlePos[particles[i].traceLength][0] = particles[i].originalPos[0];
            particles[i].particlePos[particles[i].traceLength][1] = particles[i].originalPos[1];
            particles[i].particlePos[particles[i].traceLength][2] = particles[i].originalPos[2];

            particles[i].traceLength++;
        }
    };

    AirConditioner& ac = solver->airconditioners[0];
    int x, y, z;
    x = ac.pos[0];
    y = ac.pos[1];
    z = ac.pos[2];

    if (visible == 0) {

        if (ac.acType == ACTYPE::Ceiling) {
            for (int i = 0; i <= 16; i++) {
                int ii = 17 * 0 + i;
                particles[ii].particlePos[0][0] = particles[ii].originalPos[0] = x - 4 + 0.5f;
                particles[ii].particlePos[0][1] = particles[ii].originalPos[1] = y - 0.2f;
                particles[ii].particlePos[0][2] = particles[ii].originalPos[2] = z + i / 4.0f - 2 + 0.5f;
                particles[ii].traceLength = 1;

                ii = 17 * 1 + i;
                particles[ii].particlePos[0][0] = particles[ii].originalPos[0] = x + 4 + 0.5f;
                particles[ii].particlePos[0][1] = particles[ii].originalPos[1] = y - 0.2f;
                particles[ii].particlePos[0][2] = particles[ii].originalPos[2] = z + i / 4.0f - 2 + 0.5f;
                particles[ii].traceLength = 1;

                ii = 17 * 2 + i;
                particles[ii].particlePos[0][0] = particles[ii].originalPos[0] = x + i / 4.0f - 2 + 0.5f;
                particles[ii].particlePos[0][1] = particles[ii].originalPos[1] = y - 0.2f;
                particles[ii].particlePos[0][2] = particles[ii].originalPos[2] = z - 4 + 0.5f;
                particles[ii].traceLength = 1;

                ii = 17 * 3 + i;
                particles[ii].particlePos[0][0] = particles[ii].originalPos[0] = x + i / 4.0f - 2 + 0.5f;
                particles[ii].particlePos[0][1] = particles[ii].originalPos[1] = y - 0.2f;
                particles[ii].particlePos[0][2] = particles[ii].originalPos[2] = z + 4 + 0.5f;
                particles[ii].traceLength = 1;
            }
            visible = 68;
        }
    }

    //else if (ac.acType == ACTYPE::Tower) {
    //    if (visible < numParticles - 32 && frame >= 100) {
    //        if (ac.direction % 2 == 0 && ac.direction < 8) {
    //            bool isXdir = ac.direction % 4 == 0;
    //            int flipped = isXdir ? (ac.direction == 0 ? 1 : -1) : (ac.direction == 6 ? 1 : -1);

    //            //circle
    //            {
    //                for (int i = 0; i < 12; i++) {
    //                    float angle = i / 12.f * 3.141592f / 180.f;

    //                    particles[visible].particlePos[0][0] = particles[visible].originalPos[0] = x + 0.5f;
    //                    particles[visible].particlePos[0][1] = particles[visible].originalPos[1] = 17.5f + 1.5f * cos(i) + 0.5f;
    //                    particles[visible].particlePos[0][2] = particles[visible].originalPos[2] = z - flipped * 0.5f + 1.5f * sin(i) + 0.5f;
    //                    particles[visible].traceLength = 1;
    //                    visible++;
    //                }
    //            }

    //            //guard
    //            for (int i = 5; i <= 14; i++) {


    //                particles[visible].particlePos[0][0] = particles[visible].originalPos[0] = x + 0.5f;
    //                particles[visible].particlePos[0][1] = particles[visible].originalPos[1] = i + 0.5f;
    //                particles[visible].particlePos[0][2] = particles[visible].originalPos[2] = z - 2 * flipped + 0.5f;
    //                particles[visible].traceLength = 1;
    //                visible++;

    //                particles[visible].particlePos[0][0] = particles[visible].originalPos[0] = x + 0.5f;
    //                particles[visible].particlePos[0][1] = particles[visible].originalPos[1] = i + 0.5f;
    //                particles[visible].particlePos[0][2] = particles[visible].originalPos[2] = z + 1 * flipped + 0.5f;
    //                particles[visible].traceLength = 1;
    //                visible++;

    //                //particles[visible].particlePos[0][0] = particles[visible].originalPos[0] = x + hole + 0.5f;
    //                //particles[visible].particlePos[0][1] = particles[visible].originalPos[1] = i + 0.5f;
    //                //particles[visible].particlePos[0][2] = particles[visible].originalPos[2] = z - 0.7f * flipped + 0.5f;
    //                //particles[visible].traceLength = 1;
    //                //visible++;

    //                //particles[visible].particlePos[0][0] = particles[visible].originalPos[0] = x + hole + 0.5f;
    //                //particles[visible].particlePos[0][1] = particles[visible].originalPos[1] = i + 0.5f;
    //                //particles[visible].particlePos[0][2] = particles[visible].originalPos[2] = z + 0 * flipped + 0.5f;
    //                //particles[visible].traceLength = 1;
    //                //visible++;
    //            }
    //        }

    //        //reset
    //        frame = 0;

    //    }

    //}

    glColor3f(1.f, 1.f, 1.f);
    glBegin(GL_LINES);
    //glBegin(GL_POINTS);

    for (int i = 0; i < visible; i++) {
        if (particles[i].age > lifespan) continue;

        for (int j = 1; j < particles[i].traceLength; j++) {

            positionX = particles[i].particlePos[j][0] * h - h - scale / 2.f;
            positionY = particles[i].particlePos[j][1] * h - h;
            positionZ = particles[i].particlePos[j][2] * h - h - scale / 2.f;

            prevPositionX = particles[i].particlePos[j - 1][0] * h - h - scale / 2.f;
            prevPositionY = particles[i].particlePos[j - 1][1] * h - h;
            prevPositionZ = particles[i].particlePos[j - 1][2] * h - h - scale / 2.f;

            glVertex3f(positionX, positionY, positionZ);
            glVertex3f(prevPositionX, prevPositionY, prevPositionZ);

            //glVertex3f(positionX, positionY, positionZ);
        }
    }
    glEnd();




}


void Visualizer::wuDisplay()
{
    if (!solver->isInitialized())
        return;

    solver->getGridSize(nx, ny, nz);
    maxN = nx > ny ? (nx > nz ? nx : nz) : (ny > nz ? ny : nz);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);

    glEnable(GL_ALPHA_TEST);
    glAlphaFunc(GL_GREATER, 0);

    glViewport(windowOffsetW, windowOffsetH, windowWidth, windowHeight);

    glPushMatrix();
    if (FPP)
    {
        glRotatef(rotationX, 1.0f, 0, 0);
        glRotatef(rotationY, 0, 1.0f, 0);
        glTranslatef(translationX, translationY, translationZ);
    } else {
        glLoadIdentity();
        glTranslatef(0.0f, -0.3f, -3.0f);
        glTranslatef(translationX, translationY, translationZ);
        glRotatef(rotationX, 1.0f, 0, 0);
        glRotatef(rotationY, 0, 1.0f, 0);
    }

    if (draw)
    {
        advectParticlesRK4();

        if (drawObstacle)
            DrawObstacle();

        if (drawMode % 3 == 0)
            wuDrawVelocity();
        else if (drawMode % 3 == 1)
            //wuDrawStreamsurfaceWithSplit();
            wuDrawDensity();
        else if (drawMode % 3 == 2)
            //wuDrawDensity();
            //wuDrawStreakline();
            //wuDrawStreamline(true);
            //wuDrawStreamsurfaceWithSplit();
            wuDrawStreamsurfaceWithLines();
            
    }

    if (drawAxis)
        wuDrawGrid();

    glPopMatrix();

    //glutPostRedisplay();
}

void Visualizer::wuReshape(GLint offsetW, GLint offsetH, GLint width, GLint height)
{
    windowWidth = width;
    windowHeight = height;
    windowOffsetW = offsetW;
    windowOffsetH = offsetH;

    glViewport(offsetW, offsetH, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (drawPerspective)
    {
        gluPerspective(45.0, (float)width / height, 0.001, 100.0);
    }
    else
    {
        float ratio = float(width) / height;
        if (ratio > 1)
            glOrtho(-ratio / 2.f - 1.f, 1.f + ratio / 2.f, -1, 1, 0.01, 10000);
        else
            glOrtho(-1, 1, -ratio / 2.f - 1.f, 1.f + ratio / 2.f, 0.01, 10000);
    }

    glMatrixMode(GL_MODELVIEW);
    //glLoadIdentity();
    //if(!FPP)
    //    glTranslatef(0.0f, -0.3f, -3.0f);
    //glTranslatef(translationX, translationY, translationZ);
}

void Visualizer::TemperatureDistribution(const float* t)
{
    int st[11] = { 0 };
    int N3 = int(nx * ny * nz);

    int c = 5;
    int xs = int(nx / 2 - c);
    int xf = int(nx / 2 + c);

    int zs = int(nz / 2 - c);
    int zf = int(nz / 2 + c);

    printf("////////////////////////////////////////////////////////////////////////////////////////////\n");
    printf("----------------- x = %d ~ %d  // z = %d ~ %d \n", xs, xf, zs, zf);

    for (int x = xs; x < xf; x++)
    {
        for (int y = 0; y < ny; y++)
        {
            for (int z = xs; z < xf; z++)
            {
                float tp = solver->getDensityGPU(x, y, z);
                if (tp < 0.05f)
                    continue;

                int norm = (int)(tp * 10.f);
                if (norm > 10)
                    norm = 10;
                st[norm]++;
            }
        }
    }

    int totalCnt = 0;
    for (int i = 0; i < 11; i++)
        totalCnt += st[i];

    printf("total count : %d\n", totalCnt);

    for (int i = 0; i < 11; i++)
    {
        printf("%.1f (%4.2f : %5d)\t : ", i / 10.f, float(st[i]) / totalCnt * 100.f, st[i]);

        int scaled = st[i] * 50 / totalCnt;
        for (int j = 0; j < scaled; j++)
        {
            printf("*");
        }
        printf("\n");
    }
}

void Visualizer::wuSpecialKeys(int key, int x, int y)
{
    MachineSetting& mc = setting->machineSetting[0];
    switch (key)
    {
    case GLUT_KEY_LEFT:
        mc.ventDir--;
        break;
    case GLUT_KEY_RIGHT:
        mc.ventDir++;
        break;

    case GLUT_KEY_UP:
        mc.ventSpeed++;
        break;

    case GLUT_KEY_DOWN:
        mc.ventSpeed--;
        break;
    }

    if (mc.ventSpeed < 1) mc.ventSpeed = 1;
    if (mc.ventSpeed > 10) mc.ventSpeed = 10;
    solver->setACVentMag(0, mc.ventSpeed);

    //if (mc.ventDir < 1) mc.ventDir = 1;
    //if (mc.ventDir > 4) mc.ventDir = 4;
    solver->setACVentDir(0, mc.ventDir);

    //AirConditioner ac0 = solver->airconditioners[0];
    //if (ac0.acType == ACTYPE::Ceiling) {
    //    if (mc.ventDir == 1) {
    //        solver->setACWindType(0, WINDTYPE::POWER);
    //        ac0.ventSpeed = 3;
    //    }
    //    else if (mc.ventDir == 2) {
    //        solver->setACWindType(0, WINDTYPE::FOREST);
    //        ac0.ventSpeed = 2;
    //    }
    //    else if (mc.ventDir == 3) {
    //        solver->setACWindType(0, WINDTYPE::AUTO);
    //        ac0.ventSpeed = 3;
    //    }
    //    else if (mc.ventDir == 4) {
    //        solver->setACWindType(0, WINDTYPE::AIRGUIDE);
    //        ac0.ventSpeed = 3;
    //    }
    //    else if (mc.ventDir == 5) {
    //        solver->setACWindType(0, WINDTYPE::HIGHCEILING);
    //        ac0.ventSpeed = 3;
    //    }
    //    else if (mc.ventDir == 6) {
    //        solver->setACWindType(0, WINDTYPE::FOCUS);
    //        ac0.ventSpeed = 3;
    //    }
    //}

    printf("\ndir : %d   mag : %d\n", mc.ventDir, mc.ventSpeed);
    glutPostRedisplay();
}

void Visualizer::wuKeyBoard(unsigned char key, int x, int y)
{
    //printf("%c\n", key);
    if (key == 'i' || key == 'I')
    {
        //auto name = solver->saveOcc();
        //printf("save occupancy : %s\n", name);
        //delete name;
    }
    else if (key == 'c' || key == 'C')
    {
        printf("\n -- simulation reset\n");

        pause = false;
        pauseflag = true;
        total = 0.f;
        oldTotal = 0.f;
        solver->stop();
        solver->reset(false);
        solver->start();
        logFrameStep = 0;
    }
    else if (key == 'j' || key == 'J')
    {
        printf("\n -- generate test occupnacy\n");

        solver->generateObstacle();
    }
    else if (key == 'v' || key == 'V')
    {
        drawMode++;
        visible = 0;
        if (drawMode%3==0)
            printf("\n -- draw mode : velocity\n");
        else if(drawMode % 3 == 1)
            printf("\n -- draw mode : streamsurface\n");
        else if(drawMode % 3 == 2)
            printf("\n -- draw mode : temperature\n");
    }
    else if (key == 'b' || key == 'B')
    {
        draw = !draw;
        if (draw)
            printf("\n -- draw switch : true\n");
        else
            printf("\n -- draw switch : false\n");
    }
    else if (key == ' ')
    {
        pause = !pause;
        if (pause) {
            printf("\n -- simulation pause...\n");
            //solver->stop();
        }
        else {
            printf("\n -- simulation start...\n");
            //solver->start();
        }
        pauseflag = !pause;
    }
    else if (key == 'a' || key == 'A')
    {
        rotationY -= 15;
    }
    else if (key == 'd' || key == 'D')
    {
        rotationY += 15;
    }
    else if (key == 'w' || key == 'W')
    {
        rotationX -= 10;
    }
    else if (key == 's' || key == 'S')
    {
        rotationX += 10;
    }
    else if (key == 'e' || key == 'E')
    {
        rotationX = 30;
        rotationY = -45;
    }
    else if (key == 'x' || key == 'X')
    {
        drawObstacle = !drawObstacle;
        if (drawObstacle)
            printf("\n Enable drawing obstacle\n");
        else
            printf("\n Disable drawing obstacle\n");
    }
    else if (key == 't' || key == 'T')
    {
        //drawTemperature = !drawTemperature;
        //if (drawTemperature)
        //    printf("\n Enable drawing temperature\n");
        //else
        //    printf("\n Disable drawing temperature\n");
    }
    else if (key == 'p' || key == 'P')
    {
        //auto tem = solver->getDensityGPUN();
        //TemperatureDistribution(tem);

        drawPerspective = !drawPerspective;
        if (drawPerspective)
            printf("\n Enable perspective view\n");
        else
            printf("\n Enable orthogonal view\n");

        wuReshape(windowWidth, windowHeight, windowOffsetW, windowOffsetH);
    }
    else if (key == 'o' || key == 'O')
    {
        //solver->outletSwitch();
        //if (solver->getOutletState())
        //    printf("\n -- outlet state : single\n");
        //else
        //    printf("\n -- outlet state : quadrilateral\n");
    }
    else if (key == 'l' || key == 'L')
    {
        drawAxis = !drawAxis;

        //solver->swingSwitch();
        //if (solver->getSwingState())
        //    printf("\n -- swing state : true\n");
        //else
        //    printf("\n -- swing state : false\n");
    }
    else if (key == 'k' || key == 'K')
    {
        solver->fastPlay();
    }
    else if (key == 'u' || key == 'U')
    {
        debug = !debug;
        if (debug)
        {
            solver->start();

            printf("\n -- debug mode : true\n");

            auto now = std::chrono::system_clock::now();
            auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch());
            sprintf_s(filename, "%I64d.txt", ms.count());
            OpenFile(filename);
        }
        else
        {
            solver->stop();
            //solver->reset();
            printf("\n -- debug mode : false\n");
            pause = false;
            pauseflag = true;
        }
    }
    else if (key == 'q' || key == 'Q') {
        int dir = solver->airconditioners[0].direction;
        dir++;
        if (dir >= 8)dir -= 8;
        solver->setACDirection(0,dir);
    }
    else if (key == 49) // '1'
    {
        solver->setACType(0, 1);
        //solver->setACPosition(0, int(ny) / 2, nz / 2);
        //solver->setACPosition(nx-1, int(ny) / 2, nz / 2);
        //solver->setACPosition(nx/2, int(ny) / 2, 0);
        solver->setACPosition(0, nx / 2, int(ny) / 2, nz - 1);

    }
    else if (key == 50) // '2'
    {
        solver->setACType(0, 2);
        solver->setACPosition(0, int(nx) / 2, int(ny) - 1, int(nz) / 2);
    }
    else if (key == 51) // '3' 
    {
        solver->setACType(0, 3);
        solver->setACPosition(0, 0, int(ny) * 2 / 3, int(nz) / 2);

    }
    else if (key == 52) // '4' 
    {
        drawThick--;
        if (drawThick < 5)
            drawThick = 5;;
    }
    else if (key == 53) // '5' 
    {
        drawThick++;
        if (drawThick > int(nz) - 3)
            drawThick = int(nz) - 3;
    }
    else if (key == 54) // '6' 
    {
        drawVolStart -= 1;
        if (drawVolStart < 1)
            drawVolStart = 1.f;
    }
    else if (key == 55) // '7' 
    {
        drawVolStart += 1;
        if (drawVolStart > nz - 5)
            drawVolStart = float(nz - 5);
    }
    else if (key == 56) // '8' 
    {
        if (drawScale > 0.1f)
            drawScale -= 0.1f;
        printf("\n Draw Scale : %2.1f\n", drawScale);

        //if (drawY > 0)
        //    drawY--;
        //printf("Draw Y : %d\n", drawY);
    }
    else if (key == 57) // '9' 
    {
        drawScale += 0.1f;
        printf("\n Draw Scale : %2.1f\n", drawScale);

        //if (drawY < ny-1)
        //    drawY++;
        //printf("Draw Y : %d\n", drawY);
    }
    else if (key == 'h')
    {
        help();
    }
    else if (key == 'z')
    {
        setting->envSetting.saveData = !setting->envSetting.saveData;

        if (setting->envSetting.saveData)
            printf("\n -- save simulation data : on\n");
        else
            printf("\n -- save simulation data : off\n");
    }
    else if (key == 'r')
    {
        delete solver;
        //Setting::Instance()->LoadSetting("settings/setting.json");
        setting->LoadSetting(setting->settingFileLocation.c_str());
        solver = new Solver(setting);
        total = 0.f;
        oldTotal = 0.f;
        preAvgTemp = setting->eulerianSetting.curTem;

        //if (Setting::Instance()->envSetting.saveData)
        //{
        //    debug = false;
        //}
        //else
        //{
        //    debug = true;
        //    solver->start();
        //}
    }
}

void Visualizer::setSettingVals(Setting* s)
{
    setting = s;

    minTem = s->eulerianSetting.minTem;
    maxTem = s->eulerianSetting.maxTem;
    srcTem = s->eulerianSetting.srcTem;
    curTem = s->eulerianSetting.curTem;

    nx = s->eulerianSetting.gridX;
    ny = s->eulerianSetting.gridY;
    nz = s->eulerianSetting.gridZ;
    //saveSimData = s->envSetting.saveData;

    //drawSecond = s->envSetting.drawSecond;
    //saveInterval = s->envSetting.saveInterval;
    //temHistory = s->envSetting.temHistory;
    //endTime = s->envSetting.endTime;

}

void Visualizer::wuIdle(void* a, int simNo)
{
    if (!solver->isInitialized() || !solver->isPlay())
    {
        Sleep(1);
        return;
    }

    if (debug && !pause)
    {
        maxN = nx > ny ? (nx > nz ? nx : nz) : (ny > nz ? ny : nz);
        auto& interval = setting->envSetting.saveInterval;
        if (setting->envSetting.saveData) //&& (framestep % interval == 0))
        {
            char fn[32];
            sprintf(fn, "%d.txt", solver->getFramestep());
            saveData(fn);
            AppendData(filename, solver->getSimulationTime(), setting->envSetting.scale);
        }
    }

    //solver->addSource();

    if (!pause)
    {
        //static auto tick = GetTickCount64();
        if (pauseflag) {

            solver->getGridSize(nx, ny, nz);
            maxN = nx > ny ? (nx > nz ? nx : nz) : (ny > nz ? ny : nz);

            tick = GetTickCount64();
            pauseflag = false;
        }
        if (solver->DidReset()) {
            tick = GetTickCount64();
            oldTotal = 0;
            logFrameStep = 0;
        }
        total = oldTotal + (GetTickCount64() - tick) * 0.001f;

        if (solver->getFramestep() >= logFrameStep)
        {
            logFrameStep += setting->envSetting.saveInterval;
            if (solver->getFramestep() != prevframe) {

                prevframe = solver->getFramestep();


                if (setting->envSetting.temHistory)
                    solver->pushData();

                auto timet = getTime();
                const auto avgTem = solver->getDensityAvgGPU();

                if (false) //avgTem < Setting::Instance()->machineSetting.target)
                {
                    solver->setACVentMag(0, 1);
                    for (int i = 0; i < solver->numAC; i++)
                        solver->airconditioners[i].ventSpeed = 1;
                    pause = true;
                    pauseflag = false;
                    printf("\n -- Simulation reached the target\n");
                }

                auto fps = solver->getFPS();
                int iter = (setting->envSetting.drawSecond ? int(fps) : 1);

                char str[256];
                sprintf_s(str, "Sim %d: [%s] %5dth Simtime: %4.1lfsec Total: %4.1fsec  fps: %3.1f avgT : %2.2f (%.3f)",
                    simNo, timet.c_str(), solver->getFramestep(), solver->getSimulationTime(), total, fps, avgTem, avgTem - preAvgTemp);
                //sprintf_s(str, "Sim %d: [%s] %5dth Simtime: %4.1lfsec Total: %4.1fsec  fps: %3.1f avgT : %2.2f (%.3f) locT : %2.2f",
                //    simNo, timet.c_str(), solver->getFramestep(), solver->getSimulationTime(), total, fps, avgTem, avgTem - preAvgTemp, localT);

                printf("\n%s", str);
                preAvgTemp = avgTem;

                if (simlog)
                {
                    float localT = getLocalTem(); // norm2tem(solver->getDensityGPU(px + 1, py + 1, pz + 1));

                    FILE* f;
                    char fn[512];
                    sprintf_s(fn, "%s_%d_result.txt", simlogname, ventDir);
                    if (fopen_s(&f, fn, "a") == 0)
                    {
                        char str[256];
                        fprintf(f, "%2.2f\t%2.2f\n", avgTem, localT);
                    }

                    fclose(f);

                }

                bool endSim = (setting->envSetting.endTime > 0.00001f) ? (solver->getSimulationTime() > setting->envSetting.endTime) : false;
                if (endSim) // || avgTem < (Setting::Instance()->machineSetting.tarTem))
                {
                    exit(0);
                }
            }
        }

        for (int i = 0; i < 10; i++)
        {
            //realtime
            if (realTime)
                if (solver->getSimulationTime() > total)
                    break;
            if (solver->getFramestep() >= logFrameStep)
                break;
            solver->addSource();
            solver->update();
        }

    }
    else {
        if (!pauseflag) {
            oldTotal = oldTotal + (GetTickCount64() - tick) * 0.001f;
            tick = GetTickCount64();
            pauseflag = true;
        }
    }
}

