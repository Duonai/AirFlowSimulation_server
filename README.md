# AirFlowSimulation_server

## Fluids/visualizer.cpp

서버 프로그램 상에서 기류 시뮬레이션의 형태를 가시화 하는 코드입니다.

- `void DrawVelocity()`
  - 기류 시뮬레이션 solver에서 계산한 기류의 속도, 방향 벡터를 `getVelocity`함수로 받아와서 기류의 형태를 선의 형태로 렌더링 합니다.
  - 아래 쪽으로 향하는 바람은 하얀색, 위를 향하면 초록색으로 렌더링 하여 디버깅에 용이하게 하였습니다.
 
- `void DrawObstacle()`
  - 모바일 클라이언트에서 생성해서 전송한 gird cell의 장애물 상태 여부를 시뮬레이션 solver의 getObstacle함수를 <br/> 사용해 장애물이라 확인된 grid cell들을 정육면체 오브젝트를 렌더링하여 장애물을 표현합니다.

## Fluids/TCPServer.cpp

모바일 클라이언트 기기와 통신을 위한 코드입니다.

통신 방식은 모바일 클라이언트와 비슷합니다.
