# demo
library(VCM)
# 设定样本数
n=400
# 生成z
z <- runif(n)
# 生成x
v <- rnorm(n)
x <- cumsum(v)
# theta函数
theta=function(t){(t-0.25)^2}
# 生成u
u <- rnorm(n)
# 生成y
y <- x*theta(z)+u
# 指定lambda建模
model <- val(x,z,y,p=3,nknots=20,lambda = 0.01)
# 不指定lambda建模
model <- val(x,z,y,p=3,nknots=20) #运行用时较久
# 具体参数可以查看帮助文档
help(val)

# 返回模型的lambda值
model$lambda
# 返回模型的系数
model$coef
# 画出实际变系数和其估计值散点图
plot(theta(z),model$theta)

#检验是否应该用这模型
LRTval(x, z, y, p=1,nknots=20,df=1)
help(LRTval)

#预测
z2 <- runif(100)
v2 <- rnorm(100)
x2 <- cumsum(v2)
pre <- valpredict(model,newx = x2,newz = z2)
head(pre)
