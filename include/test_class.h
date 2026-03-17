#ifndef TEST_CLASS_H_
#define TEST_CLASS_H_

class TestClass {
 public:
  TestClass();
  ~TestClass();

  void DoNothing() const;

 private:
  int value_;
};

#endif  // TEST_CLASS_H_
