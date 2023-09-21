#include "visualizer.h"
//#include "TCPServer.h"



Visualizer v[SOLVER_NUM];
Solver* solvers[SOLVER_NUM];
TCPServer server;


const int numbers[10] = { 0,1,2,3,4,5,6,7,8,9 };
const int displaygridH[12] = { 1,1,2,2,2,2,2,3,3,3,3,3 };
const int displaygridV[12] = { 1,1,1,2,2,3,3,3,3,3,4,4 };

int solverCnt = 1;
int windowWidth = 800;
int windowHeight = 800;

void Display()
{
    glClear(GL_COLOR_BUFFER_BIT);
    for (int i = 0; i < solverCnt; i++)
        v[i].wuDisplay();

    glutSwapBuffers();

    glutPostRedisplay();
}
void Reshape(GLint width, GLint height)
{
    int s = 0;

    int scale = displaygridH[solverCnt] > displaygridV[solverCnt] ? displaygridH[solverCnt] : displaygridV[solverCnt];

    int W = width / scale;
    int H = height / scale;

    int offsetH = width * (scale - displaygridH[solverCnt]) / (scale * 2);
    int offsetV = height * (scale - displaygridV[solverCnt]) / (scale * 2);

    for (int i = 0; i < displaygridV[solverCnt]; i++) {
        for (int j = 0; j < displaygridH[solverCnt]; j++) {
            v[s].wuReshape(j * W + offsetH, i * H + offsetV, W, H);

            s++;
            if (s >= solverCnt)
                break;
        }
        if (s >= solverCnt)
            break;
    }
    //if (solverCnt > 1) {
    //    v[0].wuReshape(0, 0, width / 2, height / 2);
    //    v[1].wuReshape(width / 2, 0, width / 2, height / 2);
    //    v[2].wuReshape(0, height / 2, width / 2, height / 2);
    //}
    //else {
    //    v[0].wuReshape(0, 0, width, height);
    //}
}
void SpecialKeys(int key, int x, int y)
{
    for (int i = 0; i < solverCnt; i++) {
        v[i].wuSpecialKeys(key, x, y);
    }
}
void KeyBoard(unsigned char key, int x, int y)
{
    for (int i = 0; i < solverCnt; i++) {
        v[i].wuKeyBoard(key, x, y);
    }
}
void Idle(void* a)
{
    int solvernum = *((int*)a);
    while (true)
    {
        v[solvernum].wuIdle((void*)0, solvernum);
        //v[0].wuIdle((void*)0, 0);
        //v[1].wuIdle((void*)0, 1);
    }
}

int main(int argc, char* argv[]) {

    srand(time(0));
    for (int i = 0; i < solverCnt; i++) {
        Setting* setting = new Setting();
        if (i == 0)
            setting->LoadSetting("settings/setting_manual_tower.json");
        else if (i == 1)
            setting->LoadSetting("settings/setting2.json");
        else if (i == 2)
            setting->LoadSetting("settings/setting3.json");
        else
            setting->LoadSetting("settings/setting4.json");

        v[i].setSettingVals(setting);
        v[i].preAvgTemp = setting->eulerianSetting.curTem;

        v[i].solver = new Solver(setting);
        v[i].drawScale = setting->envSetting.drawScale;
        solvers[i] = v[i].solver;
    }
    server.setSolver(solvers, solverCnt);

    if (1)
        try {
        server.open();
        thread tcpServer(&TCPServer::start, &server);
        tcpServer.detach();
    }
    catch (exception& ex) {
        printf("-- server open error\n%s\n", ex.what());
    }

    if(false)
    {
        //load BIM
        //set value
        //start simulation
        //print out log
    }



    glutInit(&argc, argv);

    for (int i = 0; i < solverCnt; i++) {
        v[i].wuInitialize();
    }



    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(windowWidth, windowHeight);
    Reshape(windowWidth, windowHeight);
    glutCreateWindow("Media Lab | Korea University");

    try {
        glutDisplayFunc(Display);
    }
    catch (exception& ex) {
        printf("-- opengl display error\n%s\n", ex.what());
    }

    glutReshapeFunc(Reshape);
    glutSpecialFunc(SpecialKeys);
    glutKeyboardFunc(KeyBoard);

    for (int i = 0; i < solverCnt; i++) {
        try {
            _beginthread(Idle, 0, (void*)&numbers[i]);
        }
        catch (exception& ex) {
            printf("-- opengl idle error\n%s\n", ex.what());
        }
    }
    glutMainLoop();
    return 0;

}