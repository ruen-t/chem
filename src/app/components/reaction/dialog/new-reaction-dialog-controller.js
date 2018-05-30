angular
    .module('app')
    .controller('DialogController', DialogController);
function DialogController($scope, $mdDialog, $http) {
    $scope.description = "";
    $scope.hide = function () {
        $mdDialog.hide();
    };

    $scope.cancel = function () {
        $mdDialog.cancel();
    };

    $scope.create = function () {
        $mdDialog.hide();
        var name = $scope.name;
        var smart = $scope.smart;
        var description = $scope.description;
        var payload = new FormData();
        var reaction = { name, smart, description };
        payload.append("reaction", JSON.stringify(reaction));
        $http({
            url: REACTION_API,
            method: 'POST',
            data: payload,
            headers: { 'Content-Type': undefined },
            transformRequest: angular.identity
        }).then(function (response) {
            console.log(response);
            swal('Successfully add new reaction', '', 'success')
        }, function (error) {
            swal('Fail to add new reaction', '', 'error')
        })

    };
}