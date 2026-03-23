#ifndef ACKROYD_RIYAIT_POINT_H_
#define ACKROYD_RIYAIT_POINT_H_

namespace ar
{
class Point
{
public:
  Point(const double x, const double y);

private:
  const double x_;
  const double y_;
};
} // namespace ar

#endif // ACKROYD_RIYAIT_POINT_H_
